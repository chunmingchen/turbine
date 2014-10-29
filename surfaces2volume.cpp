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
#include <file/nrrd.h>
#include <macros.h>

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


//run: ~/Project/turbine/build/surfaces2volume surface.dvt.list
#define EXTRACT_SURFACE

int dim[3]={150, 70, 36}; // -1 from the original sizes
const int xstart=40, xend=105; // end not included
const int ystart=10, yend=60;

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

    FILE *fout = fopen("out.raw", "wb");
    FILE *fstat = fopen("count.txt", "wt");

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
#ifdef EXTRACT_SURFACE
        vtkNew<vtkPoints> pointAry;
        vtkNew<vtkCellArray> cellAry;
#endif

        float *p = (float *)dvt_ary->GetVoidPointer(0);
        int n = data->GetNumberOfPoints();

        //vector<float> out_ary;
        vtkNew<vtkFloatArray> out_ary;
        out_ary->SetName("d_tangent_vel");
        int x,y,b;
        for (b=0; b<dim[2]; b++)
        {
            int count=0;
            for (y=ystart ;y<yend; y++)
                for (x=xstart; x<xend; x++)
                {
                    int idx = x+dim[0]*(y+dim[1]*b);
                    out_ary->InsertNextTuple1(p[idx]);
                    if (p[idx] < 0)
                        count ++;
#ifdef EXTRACT_SURFACE
                    pointAry->InsertNextPoint(data->GetPoint(idx));
                    if (x>xstart && y>ystart) {
                        vtkNew<vtkQuad> quad;
                        quad->GetPointIds()->SetId(0, xyz_id_surface(x-1,y-1,b) );
                        quad->GetPointIds()->SetId(1, xyz_id_surface(x,y-1,b) );
                        quad->GetPointIds()->SetId(2, xyz_id_surface(x,y,b) );
                        quad->GetPointIds()->SetId(3, xyz_id_surface(x-1,y,b) );
                        cellAry->InsertNextCell(quad.GetPointer());
                    }
#endif
                }

            //printf("t=%d, b=%d, count=%d\n", i, b, count);
            if (b==dim[2]-1)
                fprintf(fstat, "%d\n", count);
            else
                fprintf(fstat, "%d, ", count);
            //// separater
            //for (x=xstart; x<xend; x++)
            //    out_ary.push_back(1e+30);
        }


        //printf("n=%d\n", n);
        // fwrite(p, n, sizeof(float), fout);

        printf("i=%d, count=%d\n", i, out_ary->GetNumberOfTuples());
        fwrite(out_ary->GetVoidPointer(0), out_ary->GetNumberOfTuples(), sizeof(float), fout);
#ifdef EXTRACT_SURFACE
        vtkNew<vtkPolyData> polyout;
        polyout->SetPoints(pointAry.GetPointer());
        polyout->GetPointData()->AddArray(out_ary.GetPointer());
        polyout->SetPolys(cellAry.GetPointer());
        string filename = JCLib::strprintf("surfaces.extract.dvt.%d.vtp", i);
        vtkNew<vtkXMLPolyDataWriter> writer;
        writer->SetInputData(polyout.GetPointer());
        //writer->SetCompressorTypeToZLib();
        writer->SetFileName(filename.c_str());
        writer->Write();
#endif
    }
    write_nrrd_3d("out.nrrd", "out.raw", xend-xstart, (yend-ystart)*dim[2], files, "float");
    printf("Output: out.nrrd, out.raw\n");

    fclose(fout);
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
