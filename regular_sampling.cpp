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
#include "vtkImageCompositeDataProbeFilter.h"
#include "system/path.h"

//#define vtkImageProbeFilter vtkProbeFilter
//#define vtkImageCompositeDataProbeFilter vtkCompositeDataProbeFilter

using namespace std;

//"/data/turbine_Stg/s35_noinj_13.80_141219_turb_6201-20601"
//#define DATA_PATH "/data/turbine_Stg/zDIR.P3D.rel.6201-11001"

vtkLineWidget *lineWidget;
vtkRenderWindow *renWin;
vtkPolyData *seeds ;

vtkSmartPointer<vtkMultiPieceDataSet> load_list(char *DATA_PATH, int id)
{
  char s[1024];

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

float RES = .005;
float datamin[3] = {-0.083, -0.509, -0.509};
float datamax[3] = {0.0914, 0.509, 0.509};
struct float3 {float x,y,z;};

int main(int argc, char **argv)
{
  printf("<data_path> <id>\n");
  if (argc<=2)
    exit(1);
  char *data_path = argv[1];
  int t = atoi(argv[2]);

  // compressor
  vtkSmartPointer<vtkDataCompressor> compressor = vtkSmartPointer<vtkZLibDataCompressor>::New();

  // load data
  vtkSmartPointer<vtkMultiPieceDataSet> mb = load_list(data_path, t);

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

  vtkSmartPointer<vtkImageCompositeDataProbeFilter> resampler = vtkSmartPointer<vtkImageCompositeDataProbeFilter>::New();
  resampler->SetSourceData(mb);
  resampler->SetInputData(samplingArray);
  resampler->SpatialMatchOn();
  resampler->Update();

  //resampler->GetOutput()->PrintSelf(cout, (vtkIndent)0);

  vtkImageData *output = vtkImageData::SafeDownCast(resampler->GetOutput());

  char filename[256];
  sprintf(filename, "sampling_res%g_%d.vti", RES, t );
  printf("Save file: %s\n", filename);
  vtkNew<vtkXMLImageDataWriter> imw;
  imw->SetFileName(filename);
  imw->SetDataModeToBinary();
  imw->SetInputData(output);
  //imw->SetInputConnection(resampler.GetPointer());
  imw->SetCompressor(compressor.GetPointer());
  imw->Write();

	return 0;
}



