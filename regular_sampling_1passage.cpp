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
#include "vtkImageProbeFilter.h"

using namespace std;

vtkSmartPointer<vtkStructuredGrid> load_list(const char *gfile, const char *qfile)
{
	char s[1024];

	int i;
  vtkSmartPointer<vtkStructuredGrid> sgrid = vtkSmartPointer<vtkStructuredGrid>::New();

  printf("xyz: [%s]   q: [%s]\n", qfile, gfile);

  // Start by loading some data.
  vtkNew<vtkMultiBlockPLOT3DReader> reader;
  reader->SetXYZFileName(gfile);
  reader->SetQFileName(qfile);
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

  sgrid->DeepCopy(reader->GetOutput()->GetBlock(0));


  //extract uvel
  //<
  printf("Extracing uvel...\n");
  vtkSmartPointer<vtkFloatArray> uvel = vtkSmartPointer<vtkFloatArray>::New();
  uvel->SetName("Uvel");

  vtkDataArray* velocity_array;
  vtkPointData* PointData;
  PointData = sgrid->GetPointData();
  velocity_array = PointData->GetArray("Velocity");
  uvel->Resize(velocity_array->GetSize()/3);
  for(int p=0;p<velocity_array->GetSize()/3;p++)
  {
    double value[3];
    velocity_array->GetTuple(p,value);
    float datavalue = value[0];
    uvel->InsertTuple1(p,datavalue);
  }
  sgrid->GetPointData()->AddArray(uvel);
  printf("Done extracting uvel\n");
  //>

  return sgrid;
}

float RES = .005;
struct float3 {float x,y,z;};


int main(int argc, char **argv)
{
  printf("<gfile> <qfile> \n");
  const char *gfile = argv[1];
  const char *qfile = argv[2];

  // compressor
  vtkSmartPointer<vtkDataCompressor> compressor = vtkSmartPointer<vtkZLibDataCompressor>::New();

  // load data
  vtkSmartPointer<vtkStructuredGrid> sgrid = load_list(gfile, qfile);

  // resample
  int extent[3];
  int d,i,j,k, b;

  // set interpolator

  double *bounds = sgrid->GetBounds();
  printf("bounds: %lf %lf %lf %lf %lf %lf\n", bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);

  for (d=0; d<3; d++)
      extent[d] = (bounds[d*2+1]-bounds[d*2])/RES+1;
  printf("Extent: %d %d %d\n", extent[0], extent[1], extent[2]);


  vtkSmartPointer<vtkImageData> samplingArray = vtkSmartPointer<vtkImageData>::New();
  //samplingArray->SetNumberOfScalarComponents(1);
  //samplingArray->SetScalarType(VTK_FLOAT);
  samplingArray->SetSpacing(RES, RES, RES);
  samplingArray->SetDimensions(extent[0], extent[1], extent[2]);
  //samplingArray->SetExtent(extent);
  double origin[3];
  origin[0] = (double)bounds[0];
  origin[1] = (double)bounds[2];
  origin[2] = (double)bounds[4];
  samplingArray->SetOrigin( origin );
  //samplingArray->AllocateScalars(VTK_FLOAT,1);
  samplingArray->ComputeBounds();
  //double *bounds = samplingArray->GetBounds();
  //printf("output bounds: %f %f %f %f %f %f\n", bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);
  //bounds = vtkDataSet::SafeDownCast(mb->GetBlock(0))->GetBounds();
  //printf("input bounds: %f %f %f %f %f %f\n", bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);

  vtkSmartPointer<vtkImageProbeFilter> resampler = vtkSmartPointer<vtkImageProbeFilter>::New();
  resampler->SetSourceData(sgrid);
  resampler->SetInputData(samplingArray);
  resampler->SpatialMatchOn();
  resampler->Update();

  //resampler->GetOutput()->PrintSelf(cout, (vtkIndent)0);

  vtkImageData *output = vtkImageData::SafeDownCast(resampler->GetOutput());

  char filename[256];
  sprintf(filename, "%s.vti", qfile, RES);
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



