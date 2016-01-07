// Reference: VTK/Examples/VisualizationAlgorithms/Python/StreamlinesWithLineWidget.py
#include <stdio.h>
#include <stdlib.h>

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
#include "vtkStructuredGridOutlineFilter.h"
#include "vtkNew.h"
#include "vtkSmartPointer.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkXMLMultiBlockDataWriter.h"
#include "vtkZLibDataCompressor.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"

//#include "system/path.h"

using namespace std;

#define DATA_PATH "./"    //"/data/flow2/Stg37/"

const int VARS=4;
char *var_list[VARS] = {"Density", "Pressure", "Temperature", "VelocityMagnitude"};

// load from file list
// omits writing the file_list file
extern vtkSmartPointer<vtkMultiBlockDataSet> load_list(vector<string> &files)
{
  string dPath = "";// Dataset::getCurrentDataDir();
  vtkSmartPointer<vtkMultiBlockDataSet> mb = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  int blocks = files.size()/2;
#pragma omp parallel for
  for (int i=0; i<blocks; i++)
  {
      // Start by loading some data.
      vtkSmartPointer<vtkMultiBlockPLOT3DReader> pl3dReader = vtkSmartPointer<vtkMultiBlockPLOT3DReader>::New();
      pl3dReader->SetXYZFileName((dPath+files[i*2]).c_str());
      pl3dReader->SetQFileName((dPath+files[i*2+1]).c_str());
      pl3dReader->SetAutoDetectFormat(1);

      pl3dReader->SetScalarFunctionNumber(100);
      pl3dReader->SetVectorFunctionNumber(200);

      char scalarVarList[] = {1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0};

      if(scalarVarList[0])
      pl3dReader->AddFunction(100); //Density
      if(scalarVarList[1])
      pl3dReader->AddFunction(110); //Pressure
      if(scalarVarList[2])
      pl3dReader->AddFunction(120); //Temperature
      if(scalarVarList[3])
      pl3dReader->AddFunction(130); //Enthalpy
      if(scalarVarList[4])
      pl3dReader->AddFunction(140); //Internal Energy
      if(scalarVarList[5])
      pl3dReader->AddFunction(144); //Kinetic Energy
      if(scalarVarList[6])
      pl3dReader->AddFunction(153); //Velocity Magnitude
      if(scalarVarList[7])
      pl3dReader->AddFunction(163); //Stagnation Energy
      if(scalarVarList[8])
      pl3dReader->AddFunction(170); //Entropy
      if(scalarVarList[9])
      pl3dReader->AddFunction(184); //Swirl

      //Vector Fields
      if(scalarVarList[10])
      pl3dReader->AddFunction(200); //Velocity
      if(scalarVarList[11])
      pl3dReader->AddFunction(201); //Vorticity
      if(scalarVarList[12])
      pl3dReader->AddFunction(202); //Momentum
      if(scalarVarList[13])
      pl3dReader->AddFunction(210); //Pressure Gradient

      pl3dReader->Update();
      vtkDataSet *current_data = vtkDataSet::SafeDownCast(pl3dReader->GetOutput()->GetBlock(0));
#pragma omp critical
      mb->SetBlock(i, current_data);

      printf("Block %d loaded\n", i);
  }
  //mb->GlobalReleaseDataFlagOn(); //!!! 11/11/14

  return mb;
}

vector<string> create_file_list(char *path, int numStep, int blockId, int numBlocks)
{
    vector<string> files;
    string init_name = "s35_noinj.r2b";
    string middle_name = ".p3d.";
    stringstream a;
    string name1,name2;
    int actualTimeStep= numStep; // replace time step 52 with spinner variable
    a<<actualTimeStep;

    //Start Writing list file
    for(int i=blockId+1;i<=blockId+numBlocks;i++) //numBlocksToLoad
    {
        stringstream b;
        b<<i;
        name1 = string(path) + "/" + init_name + b.str() + middle_name + "g" + a.str();
        name2 = string(path) + "/" + init_name + b.str() + middle_name + "q" + a.str();
        files.push_back(name1);
        files.push_back(name2);
    }
    return files;
}

int main(int argc, char **argv)
{
  printf("convertUnstructured <data_path> <Timestep>\n");
  if (argc<3)
    return 0;

  char *path = argv[1];
  int t = atoi(argv[2]);

  vector<string> files = create_file_list(path, t, 0, 36);
  vtkSmartPointer<vtkMultiBlockDataSet> mb = load_list(files);

  for (int b=0; b<36; b++)
  {
    printf("Output block %d...\n", b+1);
    vtkSmartPointer<vtkStructuredGrid> sgrid = vtkStructuredGrid::SafeDownCast( mb->GetBlock(b) );
    int i,j,k;
    FILE *fp;
    char s[1024];

    // output grids:
    sprintf(s, "t%d.r2.b%d.grid", t, b+1);
    fp = fopen(s, "wb");
    vtkPoints *points = sgrid->GetPoints();
    int len = points->GetNumberOfPoints();
    fwrite(&len, 4, 1, fp);
    for (i=0; i<points->GetNumberOfPoints(); i++)
    {
      double x[3];
      points->GetPoint(i, x);
      float f[3];
      f[0] = x[0]; f[1] = x[1]; f[2] = x[2];
      fwrite(f, sizeof(float), 3, fp);
    }
    fclose(fp);

    // output scalars
    vtkPointData *data = sgrid->GetPointData();
    for (int v = 0; v<VARS; v++)
    {
      vtkDataArray *array = data->GetArray(var_list[v]);
      sprintf(s, "t%d.r2.b%d.%s.raw", t, b+1, array->GetName());
      fp = fopen(s, "wb");
      int len = array->GetNumberOfTuples();
      fwrite(&len, sizeof(int), 1, fp);
      double range[2];
      array->GetRange(range);
      float min = range[0],
            max = range[1];
      printf("%s %f %f\n", var_list[v], min, max);
      fwrite(&min, sizeof(float), 1, fp);
      fwrite(&max, sizeof(float), 1, fp);
      fwrite(array->GetVoidPointer(0), sizeof(float), len, fp);
      fclose(fp);
    }

    if (b==0) {
      // output cells
      sprintf(s, "t%d.r2.cell", t);
      int *dim = sgrid->GetDimensions();

      fp = fopen(s, "wb");
      int len = (dim[0]-1)*(dim[1]-1)*(dim[2]-1);
      fwrite(&len, sizeof(int), 1, fp);
      for (k=0; k<dim[2]-1; k++)
        for (j=0; j<dim[1]-1; j++)
          for (i=0; i<dim[0]-1; i++)
          {
  #define IJK_ID(i,j,k) (i+dim[0]*(j+dim[1]*(k)))
            int id[8] = {
              IJK_ID(i  ,j  ,k  ),
              IJK_ID(i+1,j  ,k  ),
              IJK_ID(i,j+1  ,k  ),
              IJK_ID(i+1,j+1,k  ),
              IJK_ID(i  ,j  ,k+1),
              IJK_ID(i+1,j  ,k+1),
              IJK_ID(i  ,j+1,k+1),
              IJK_ID(i+1,j+1,k+1)
            };
            fwrite(id, sizeof(int), 8, fp);
          }
      fclose(fp);
    }
  }

  return 0;
}



