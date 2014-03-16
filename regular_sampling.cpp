// Reference: VTK/Examples/VisualizationAlgorithms/Python/StreamlinesWithLineWidget.py
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include <list>
#include <iterator>

//#include "vtkOSUFlow.h"
#include "vtkDataSet.h"
#include "vtkStructuredGrid.h"
#include "vtkMultiBlockPLOT3DReader.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkSmartPointer.h"
#include "vtkLineSource.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkStructuredGridOutlineFilter.h"
#include "vtkProperty.h"
#include "vtkLineWidget.h"
#include "vtkCommand.h"
#include "vtkCallbackCommand.h"
#include "vtkNew.h"
#include "vtkSmartPointer.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkXMLMultiBlockDataWriter.h"
#include "vtkZLibDataCompressor.h"
// interpolate
#include "vtkInterpolatedVelocityField.h"

using namespace std;

#define DATA_PATH "./" //"/data/flow2/Stg37/"

vtkLineWidget *lineWidget;
vtkRenderWindow *renWin;
vtkPolyData *seeds ;


vtkSmartPointer<vtkMultiBlockDataSet> load_list(char *list_fname)
{
	//char list_filename[1024];
	//sprintf(list_filename, "%s/list", DATA_PATH);
	FILE *fp = fopen(list_fname, "rt");
	char s[1024];
	fgets(s, 1024, fp);

	int blocks = atoi(s);
	int i;
	vtkSmartPointer<vtkMultiBlockDataSet> mb = vtkMultiBlockDataSet::New();

	for (i=0; i<blocks; i++)
	{
		char file1[1024], file2[1024]; // q, xyz
		fgets(s, 1024, fp);
		*strchr(s, '\n')=0; // remove the last new-line
		sprintf(file1, "%s%s", DATA_PATH, s);
		fgets(s, 1024, fp);
		*strchr(s, '\n')=0; // remove the last new-line
		sprintf(file2, "%s%s", DATA_PATH, s);
		printf("xyz: [%s]   q: [%s]\n", file1, file2);

		// Start by loading some data.
		vtkNew<vtkMultiBlockPLOT3DReader> reader;
		reader->SetXYZFileName(file1);
		reader->SetQFileName(file2);
		reader->SetScalarFunctionNumber(100);
		reader->SetVectorFunctionNumber(200);
		reader->SetAutoDetectFormat(1);

	    reader->AddFunction(100); //density
	    //reader->AddFunction(110); //pressure
	    //reader->AddFunction(120); //temp
	    //reader->AddFunction(130); //enthalpy
	    //reader->AddFunction(140); //internal energy
	    //reader->AddFunction(144); //kinetic energy
	    //reader->AddFunction(153); //vel magnitude
	    //reader->AddFunction(163); //stagnation energy
	    //reader->AddFunction(170); //entropy
	    //reader->AddFunction(184); //swirl
	    //reader->AddFunction(211); //vorticity magnitude

	    //available vector fields in the data
	    reader->AddFunction(200); //velocity
	    //reader->AddFunction(201); //vorticity
	    //reader->AddFunction(202); //momentum
	    //reader->AddFunction(210); //pressure gradient
	    //reader->AddFunction(212); //starin rate

	    reader->Update();
		vtkDataSet *current_data = vtkDataSet::SafeDownCast(reader->GetOutput()->GetBlock(0));



		mb->SetBlock(i, current_data);

	}
	return mb;
}

float RES = .01;
float datamin[3] = {-0.083, -0.509, -0.509};
float datamax[3] = {0.0914, 0.509, 0.509};
struct float3 {float x,y,z;};

int main(int argc, char **argv)
{
	printf("Usage: convertVTK <list file> <output>\n");

    // compressor
    vtkSmartPointer<vtkDataCompressor> compressor = vtkZLibDataCompressor::New();

	// load data
	vtkSmartPointer<vtkMultiBlockDataSet> mb = load_list(argv[1]);

	// resample
	int extent[3];
	int d,i,j,k, b;

	// set interpolator

	vtkInterpolatedVelocityField *vecInterp = vtkInterpolatedVelocityField::New();
	for (b=0; b<mb->GetNumberOfBlocks(); b++)
	{
		vtkDataSet *dataset = vtkDataSet::SafeDownCast( mb->GetBlock(b) );
		double *bounds = dataset->GetBounds();
		//printf("bounds: %f %f %f %f %f %f\n", bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5], bounds[6]);

		// vector
		vecInterp->AddDataSet(dataset);
	}

	for (d=0; d<3; d++)
		extent[d] = (datamax[d]-datamin[d])/RES+1;
	printf("Extent: %d %d %d\n", extent[0], extent[1], extent[2]);
	vector<float3> vec_ary(extent[0]*extent[1]*extent[2]);

	for (k=0; k<extent[2]; k++)
	{
		for (j=0; j<extent[1]; j++)
			for (i=0; i<extent[0]; i++)
			{
				int b, res=0;
				double  coords[4];

				coords[0] = i*RES+datamin[0];
				coords[1] = j*RES+datamin[1];
				coords[2] = k*RES+datamin[2];
				coords[3] = 0;
				// query data
				double vec[3];
				res = vecInterp->FunctionValues(coords, vec) ; // returns nonzero if success

				float3 output;
				if (!res) {
					output.x = output.y = output.z = 0; //1e+10;
					//printf("[%lf %lf %lf] output: NaN\n", coords[0], coords[1], coords[2]);
				} else {
					output.x = (float)vec[0];
					output.y = (float)vec[1];
					output.z = (float)vec[2];
					//printf("output: %f %f %f\n", output.x, output.y, output.z);
				}
				vec_ary[i+extent[0]*(j+extent[1]*k)] = output;

			}
		printf("%d/%d\n", k, extent[2]);
	}
	char filename[256];
	char *namebase ;
	if (argc==2)
		namebase = argv[1];
	else
		namebase = argv[2];

	printf("output: %s.raw, %s.nhdr\n", namebase, namebase);

	sprintf(filename, "%s.raw", namebase);
	FILE *fp = fopen(filename, "wb");
	fwrite(&vec_ary[0], sizeof(float3), vec_ary.size(), fp);
	fclose(fp);

	sprintf(filename, "%s.nhdr", namebase);
	fp = fopen(filename, "wt");
	fprintf(fp, "NRRD0004\n"
			"type: float\n"
			"dimension: 3\n"
			"sizes: %d %d %d\n"
			"encoding: raw\n"
			"data file: %s.raw\n"
			"space origin: (%f,%f,%f)\n"
			"space directions: (%f,0,0) (0,%f,0) (0,0,%f)\n"
			, extent[0], extent[1], extent[2],
			namebase,
			datamin[0], datamin[1], datamin[2],
			RES, RES, RES);
	fclose(fp);

#if 0
	vtkNew<vtkXMLMultiBlockDataWriter> mbw;
	mbw->SetFileName("merged.vtm");
	mbw->SetDataModeToBinary();
	mbw->SetInputData(mb.GetPointer());
    mbw->SetCompressor(compressor);
	mbw->Write();
#endif

	return 0;
}



