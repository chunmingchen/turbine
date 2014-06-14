//C++ standard Header Files
#include <iostream>
#include <cmath>
#include <ctime>
#include <vector>
#include <sstream>
#include <string>
#include <cstring>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <algorithm>
#include <map>
#include <assert.h>
#include <system/path.h>

// VTK#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkVersion.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataReader.h>
#include <vtkPoints.h>
#include <vtkMultiPieceDataSet.h>
#include <vtkNew.h>
#include <vtkMultiBlockPLOT3DReader.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkQuad.h>

#include "VectorMatrix.h"


int dim[3]={150, 70, 36}; // -1 from the original sizes
const int xstart=40, xend=120;

void run(char *list_fname)
{
    //char list_filename[1024];
	//sprintf(list_filename, "%s/list", DATA_PATH);
	FILE *fp = fopen(list_fname, "rt");
	char s[1024];
	fgets(s, 1024, fp);

    int files = atoi(s);
    int i;

    FILE *fout = fopen("out.raw", "wb");

    for (i=0; i<files; i++)
    {
        fgets(s, 1024, fp);
        char fname[1024];
        sscanf(s, "%s", fname);

        // load with vtk
        vtkNew<vtkXMLPolyDataReader> reader;
        reader->SetFileName(fname);
        reader->UpdateInformation();
        reader->Update();

        // get data
        vtkSmartPointer<vtkPolyData> data = reader->GetOutput() ;

        // output
        vtkFloatArray* dvt_ary = vtkFloatArray::SafeDownCast(
                    data->GetPointData()->GetArray("d_tangent_vel")
                    );
        if (dvt_ary==NULL) {
            printf("cannot get d_tangent_vel\n");
            exit(1);
        }

        float *p = (float *)dvt_ary->GetVoidPointer(0);
        int n = data->GetNumberOfPoints();

        vector<float> out_ary;
        int x,y,z;
        for (z=0; z<dim[2]; z++)
            for (y=0 ;y<dim[1]; y++)
                for (x=xstart; x<xend; x++)
                {
                    int idx = x+dim[0]*(y+dim[1]*z);
                    out_ary.push_back(p[idx]);
                }


        //printf("n=%d\n", n);
        // fwrite(p, n, sizeof(float), fout);

        printf("i=%d, count=%d\n", i, out_ary.size());
        fwrite(&out_ary[0], out_ary.size(), sizeof(float), fout);

    }

    fclose(fout);
    fclose(fp);

}

int main(int argc, char **argv)
{
    printf("usage: filename \n");
    char *filename = argv[1];

	// load data
    run(filename);


}
