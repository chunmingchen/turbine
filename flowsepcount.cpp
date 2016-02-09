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
#include <jclib/system/path.h>
#include <jclib/file/nrrd.h>
#include <jclib/macros.h>

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
#include <vtkPoints.h>
#include <vtkXMLPolyDataWriter.h>

#include "VectorMatrix.h"


//run: ~/Project/turbine/build/flowsepcount surface.dvt.list
#define EXTRACT_SURFACE

int dim[3]={150, 70, 36}; // -1 from the original sizes
const int xstart=40, xend=105; // end not included
const int ystart=30, yend=61; // hub to tip (end not included)

inline int xyz_id_surface(int x, int y, int b)
{
    return (x-xstart)+(xend-xstart)*(y-ystart + (yend-ystart)*b);
}

void run(char *list_fname)
{
    //char list_filename[1024];
	//sprintf(list_filename, "%s/list", DATA_PATH);
	FILE *fp = fopen(list_fname, "rt");
	char s[1024];
	fgets(s, 1024, fp);

    int files = atoi(s);
    int i;

    FILE *fstat = fopen("flowsepcount.txt", "wt");

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

        int x,y,b;
        for (b=0; b<dim[2]; b++)
        {
            vector<float> out_ary(xend-xstart, 0.f);
            int count=0;
            for (y=ystart ;y<yend; y++)
                for (x=xstart; x<xend; x++)
                {
                    int idx = x+dim[0]*(y+dim[1]*b);
                    if (p[idx] < 0)
                        out_ary[x-xstart]  = max(out_ary[x-xstart], -p[idx]);
                }

            for (x=0; x<xend-xstart; x++)
            {
                fprintf(fstat, "%g ", out_ary[x]);
            }
            fprintf(fstat, "\n");
        }

        printf("%d\n", i);
        fflush(fstat);
    }
    fclose(fp);
    fclose(fstat);

}

int main(int argc, char **argv)
{
    printf("usage: filename \n");
    char *filename = argv[1];

	// load data
    run(filename);


}
