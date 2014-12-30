// Reference: VTK/Examples/VisualizationAlgorithms/Python/StreamlinesWithLineWidget.py
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>

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
#include "vtkMultiPieceDataSet.h"
#include "vtkXMLMultiBlockDataWriter.h"
#include "vtkZLibDataCompressor.h"
// interpolate
#include "vtkInterpolatedVelocityField.h"
#include "vtkImageData.h"
#include "vtkProbeFilter.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkCompositeDataProbeFilter.h"
#include "vtkXMLImageDataWriter.h"

using namespace std;

#define DATA_PATH DATA_DIR //"/data/turbine_Stg/s35_noinj_13.80_141219_turb_6201-20601"
//#define DATA_PATH "/data/turbine_Stg/zDIR.P3D.rel.6201-11001"

vtkLineWidget *lineWidget;
vtkRenderWindow *renWin;
vtkPolyData *seeds ;

vtkSmartPointer<vtkMultiPieceDataSet> load_list(int t)
{
	//char list_filename[1024];
	//sprintf(list_filename, "%s/list", DATA_PATH);
    //FILE *fp = fopen(list_fname, "rt");
	char s[1024];
    //fgets(s, 1024, fp);

    int blocks = 36;
	int i;
    vtkSmartPointer<vtkMultiPieceDataSet> mb = vtkSmartPointer<vtkMultiPieceDataSet>::New();

	for (i=0; i<blocks; i++)
	{
		char file1[1024], file2[1024]; // q, xyz
#if 0
		fgets(s, 1024, fp);
		*strchr(s, '\n')=0; // remove the last new-line
		sprintf(file1, "%s/%s", DATA_PATH, s);
		fgets(s, 1024, fp);
		*strchr(s, '\n')=0; // remove the last new-line
		sprintf(file2, "%s/%s", DATA_PATH, s);
#else
        int id = t*25+6201;
        sprintf(file1, "%s/s35_noinj.r2b%d.p3d.g%d", DATA_PATH, i+1, id);
        sprintf(file2, "%s/s35_noinj.r2b%d.p3d.q%d", DATA_PATH, i+1, id);
#endif
		printf("xyz: [%s]   q: [%s]\n", file1, file2);

		// Start by loading some data.
		vtkNew<vtkMultiBlockPLOT3DReader> reader;
		reader->SetXYZFileName(file1);
		reader->SetQFileName(file2);
		reader->SetScalarFunctionNumber(100);
		reader->SetVectorFunctionNumber(200);
		reader->SetAutoDetectFormat(1);

	    reader->AddFunction(100); //density
	    reader->AddFunction(110); //pressure
        //reader->AddFunction(120); //temp
	    //reader->AddFunction(130); //enthalpy
	    //reader->AddFunction(140); //internal energy
	    //reader->AddFunction(144); //kinetic energy
	    //reader->AddFunction(153); //vel magnitude
	    //reader->AddFunction(163); //stagnation energy
	    reader->AddFunction(170); //entropy
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

        //cout << current_data->GetNumberOfPoints() << endl;

		//extract uvel
		//<
        printf("Extracing uvel...\n");
		vtkSmartPointer<vtkFloatArray> uvel = vtkSmartPointer<vtkFloatArray>::New();
		uvel->SetName("Uvel");

		vtkDataArray* velocity_array;
		vtkPointData* PointData;
		PointData = current_data->GetPointData();
		velocity_array = PointData->GetArray("Velocity");
		uvel->Resize(velocity_array->GetSize()/3);
		for(int p=0;p<velocity_array->GetSize()/3;p++)
			{
				double value[3];
			 	velocity_array->GetTuple(p,value);
			 	float datavalue = value[0];
			 	uvel->InsertTuple1(p,datavalue);
			}
	    current_data->GetPointData()->AddArray(uvel);
        printf("Done extracting uvel\n");
	    //>

	    mb->SetPiece(i, current_data);


	}
	return mb;
}

float RES = .01; //.005;
float datamin[3] = {-0.083, -0.509, -0.509};
float datamax[3] = {0.0914, 0.509, 0.509};
struct float3 {float x,y,z;};

void run(int t)
{
    // compressor
    vtkSmartPointer<vtkDataCompressor> compressor = vtkSmartPointer<vtkZLibDataCompressor>::New();

    // load data
    vtkSmartPointer<vtkMultiPieceDataSet> mb = load_list(t);

    // resample
    int extent[3];
    int d,i,j,k, b;

    // set interpolator


    for (d=0; d<3; d++)
        extent[d] = (datamax[d]-datamin[d])/RES+1;
    printf("Extent: %d %d %d\n", extent[0], extent[1], extent[2]);


    vtkSmartPointer<vtkImageData> samplingArray = vtkSmartPointer<vtkImageData>::New();
    //samplingArray->SetNumberOfScalarComponents(1);
    //samplingArray->SetScalarType(VTK_FLOAT);
    samplingArray->SetSpacing(RES, RES, RES);
    samplingArray->SetDimensions(extent[0], extent[1], extent[2]);
    double origin[3];
    origin[0] = (double)datamin[0];
    origin[1] = (double)datamin[1];
    origin[2] = (double)datamin[2];
    samplingArray->SetOrigin( origin );
    samplingArray->AllocateScalars(VTK_FLOAT,1);
    samplingArray->ComputeBounds();
    double *bounds = samplingArray->GetBounds();
    printf("output bounds: %f %f %f %f %f %f\n", bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);
    //bounds = vtkDataSet::SafeDownCast(mb->GetBlock(0))->GetBounds();
    //printf("input bounds: %f %f %f %f %f %f\n", bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);

    vtkSmartPointer<vtkCompositeDataProbeFilter> resampler = vtkSmartPointer<vtkCompositeDataProbeFilter>::New();
    resampler->SetSourceData(mb);
    resampler->SetInputData(samplingArray);
    resampler->SpatialMatchOn();
    resampler->Update();

    //resampler->GetOutput()->PrintSelf(cout, (vtkIndent)0);

    vtkImageData *output = vtkImageData::SafeDownCast(resampler->GetOutput());

    char filename[256];
    sprintf(filename, "regular_r%g_%d.vti", RES, t);
    vtkNew<vtkXMLImageDataWriter> imw;
    imw->SetFileName(filename);
    imw->SetDataModeToBinary();
    imw->SetInputData(output);
    //imw->SetInputConnection(resampler.GetPointer());
    imw->SetCompressor(compressor.GetPointer());
    imw->Write();

}


int main(int argc, char **argv)
{
#pragma omp parallel for
    for (int i=0; i<576; i++)
    {
        run (i);
    }
#if 0
    printf("Usage: convertVTK <timestep_0base>\n");

    // compressor
    vtkSmartPointer<vtkDataCompressor> compressor = vtkSmartPointer<vtkZLibDataCompressor>::New();

	// load data
    vtkSmartPointer<vtkMultiPieceDataSet> mb = load_list(atoi(argv[1]));

	// resample
	int extent[3];
	int d,i,j,k, b;

	// set interpolator


	for (d=0; d<3; d++)
		extent[d] = (datamax[d]-datamin[d])/RES+1;
	printf("Extent: %d %d %d\n", extent[0], extent[1], extent[2]);

#if 1

    vtkSmartPointer<vtkImageData> samplingArray = vtkSmartPointer<vtkImageData>::New();
	//samplingArray->SetNumberOfScalarComponents(1);
	//samplingArray->SetScalarType(VTK_FLOAT);
	samplingArray->SetSpacing(RES, RES, RES);
	samplingArray->SetDimensions(extent[0], extent[1], extent[2]);
	double origin[3];
	origin[0] = (double)datamin[0];
	origin[1] = (double)datamin[1];
	origin[2] = (double)datamin[2];
	samplingArray->SetOrigin( origin );
	samplingArray->AllocateScalars(VTK_FLOAT,1);
	samplingArray->ComputeBounds();
	double *bounds = samplingArray->GetBounds();
	printf("output bounds: %f %f %f %f %f %f\n", bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);
	//bounds = vtkDataSet::SafeDownCast(mb->GetBlock(0))->GetBounds();
	//printf("input bounds: %f %f %f %f %f %f\n", bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);

    vtkSmartPointer<vtkCompositeDataProbeFilter> resampler = vtkSmartPointer<vtkCompositeDataProbeFilter>::New();
	resampler->SetSourceData(mb);
	resampler->SetInputData(samplingArray);
	resampler->SpatialMatchOn();
	resampler->Update();

    resampler->GetOutput()->PrintSelf(cout, (vtkIndent)0);

    vtkImageData *output = vtkImageData::SafeDownCast(resampler->GetOutput());

	char filename[256];
    sprintf(filename, "regular_%s_r%g.vti", argv[1], RES);
	vtkNew<vtkXMLImageDataWriter> imw;
	imw->SetFileName(filename);
	imw->SetDataModeToBinary();
	imw->SetInputData(output);
	imw->SetCompressor(compressor);
    imw->Write();
#if 0
	//vtkImageData *image = vtkImageData::New();
	//image->DeepCopy(resampler->GetImageDataOutput());
	vtkDataArray* data = resampler->GetOutput()->GetPointData()->GetScalars();
	vtkSmartPointer<vtkFloatArray> ary = vtkFloatArray::SafeDownCast(data);
	//for (int i=0; i<ary->GetNumberOfTuples(); i++) {
	//	if (ary->GetValue(i) != 0 )
	//		printf("%f, ", ary->GetValue(i));
	//}

	char filename[256];
	char *namebase ;
	if (argc==2)
		namebase = argv[1];
	else
		namebase = argv[2];

	printf("output: %s_v.raw, %s_v.nhdr\n", namebase, namebase);

	sprintf(filename, "%s.raw", namebase);
	FILE *fp = fopen(filename, "wb");
	fwrite(ary->GetVoidPointer(0), sizeof(float), ary->GetNumberOfTuples(), fp);
	fclose(fp);

	sprintf(filename, "%s_v.nhdr", namebase);
	fp = fopen(filename, "wt");
	fprintf(fp, "NRRD0004\n"
			"type: float\n"
			"dimension: 3\n"
			"sizes: %d %d %d\n"
			"encoding: raw\n"
			"data file: %s_v.raw\n"
			"space origin: (%f,%f,%f)\n"
			"space directions: (%f,0,0) (0,%f,0) (0,0,%f)\n"
			, extent[0], extent[1], extent[2],
			namebase,
			datamin[0], datamin[1], datamin[2],
			RES, RES, RES);
	fclose(fp);
#endif

#else
	vtkInterpolatedVelocityField *vecInterp = vtkInterpolatedVelocityField::New();
	for (b=0; b<mb->GetNumberOfPieces(); b++)
	{
		vtkDataSet *dataset = vtkDataSet::SafeDownCast( mb->GetPiece(b) );
		double *bounds = dataset->GetBounds();
		//printf("bounds: %f %f %f %f %f %f\n", bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5], bounds[6]);

		// vector
		vecInterp->AddDataSet(dataset);
	}
	vector<float3> vec_ary(extent[0]*extent[1]*extent[2]);

	// extract velocity
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

#endif
#if 0
	vtkNew<vtkXMLMultiBlockDataWriter> mbw;
	mbw->SetFileName("merged.vtm");
	mbw->SetDataModeToBinary();
	mbw->SetInputData(mb.GetPointer());
    mbw->SetCompressor(compressor);
	mbw->Write();
#endif
#endif

	return 0;
}



