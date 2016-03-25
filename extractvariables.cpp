// extract variables to individual files from vtk files

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <vtkNew.h>
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
#include "vtkXMLImageDataReader.h"
#include "vtkImageData.h"
#include "vtkPointData.h"
#include "jclib/system/path.h"
#include "jclib/system/cmd_arg_reader.h"

#define vsp_new(type, x) vtkSmartPointer<type> x = vtkSmartPointer<type>::New()

const int VARS=8;
const char *var_list[] = {"Pressure", "Entropy", "Density", "Lambda2", "Temperature", "Uvel", "VelocityMagnitude", "vtkValidPointMask"}; //, "Entropy"};
//const int VARS=1;
//const char *var_list[] = {"emd"};


int main(int argc, const char **argv)
{
  CmdArgReader::init(argc, argv, "-slice=-1");

  printf("Usage: <vti filename> -slice=<x slice idx>\n");
  vsp_new(vtkXMLImageDataReader, reader);
  reader->SetFileName(argv[1]);
  reader->Update();
  vtkImageData *image = vtkImageData::SafeDownCast( reader->GetOutput() );
  vtkPointData *data = image->GetPointData();

  // get ext
  const int *dim = image->GetDimensions();
  printf("Extent: %d %d %d\n", dim[0], dim[1], dim[2]);

  int slice = GET_ARG_INT("slice");

#if 0 // all variables
  for (int i=0; i<data->GetNumberOfArrays(); i++)
  {
    vtkDataArray *array = data->GetArray(i);
#else // select variables
  for (int v=0; v<VARS; v++)
  {
    vtkDataArray *array = data->GetArray(var_list[v]);

#endif
    char fname[1024];
    if (slice==-1)
      sprintf(fname, "%s_%s.raw", getFilename(argv[1]).c_str(), array->GetName());
    else
      sprintf(fname, "%s_%s_slice%d.raw", getFilename(argv[1]).c_str(), array->GetName(), slice);

    FILE *fp = fopen(fname, "wb");

    printf("Writing %s\n", fname);

    int ii,jj,kk;
    int istart = (slice>=0)?slice:0;
    int iend = (slice>=0)? slice+1: dim[0];
    for (kk=0; kk<dim[2]; kk++)
      for (jj=0; jj<dim[1]; jj++)
        for (ii=istart; ii<iend; ii++)
        {
          int idx = ii+dim[0]*(jj+dim[1]*kk);
          float x = (float)array->GetTuple1(idx);
          fwrite(&x, sizeof(float), 1, fp);
        }
    fclose(fp);
  }
  printf("Done\n");
  return 0;
}
