#include <sstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <vtkXMLMultiBlockDataReader.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkInterpolatedVelocityField.h>
#include <vtkNew.h>
#include <vtkMultiBlockPLOT3DReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkDataSet.h>
#include <vtkStructuredData.h>
#include <vtkStructuredGrid.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkXMLStructuredGridWriter.h>

int main(int argc, char ** argv)
{
    const char *label_file = argv[1];
    const char *grid_file = argv[2];

    vtkNew<vtkXMLStructuredGridReader> reader;
    reader->SetFileName(grid_file);
    reader->Update();
    vtkStructuredGrid *data = vtkStructuredGrid::SafeDownCast(reader->GetOutput());
    int *ext = data->GetExtent();
    int dim[3] = {ext[1]-ext[0]+1, ext[3]-ext[2]+1, ext[5]-ext[4]+1};
    printf("data dim: %d %d %d\n", dim[0], dim[1], dim[2]);

    int in_ext[6];
    FILE *fp = fopen(label_file, "rt");
    if (!fp) {
        perror(label_file);
        return 1;
    }
    fscanf(fp, "%d %d %d %d %d %d", &in_ext[1], &in_ext[0], &in_ext[3], &in_ext[2], &in_ext[5], &in_ext[4]);

    int i,j,k, count=0;
    vtkNew<vtkPoints> points;
    vtkNew<vtkIntArray> field;
    for (k=in_ext[4]; k<in_ext[5]; k++)
        for (j=in_ext[2]; j<in_ext[3]; j++)
            for (i=in_ext[0]; i<in_ext[1]; i++)
            {
                int idx = i+dim[0]*(j+dim[1]*k);
                double *p = data->GetPoint(idx);
                points->InsertNextPoint(p);

                int label;
                int r = fscanf(fp, "%d", &label);
                if (r!=1) {
                    perror(label_file);
                    return 1;
                }
                field->InsertNextValue(label);
            }
    field->SetName("Label");


    vtkNew<vtkStructuredGrid> outData;
    outData->SetDimensions(in_ext[1]-in_ext[0], in_ext[3]-in_ext[2], in_ext[5]-in_ext[4]);
    outData->SetPoints(points.GetPointer());
    outData->GetPointData()->AddArray(field.GetPointer());


    std::stringstream ss;
    ss << label_file << ".vts";
    printf("Output: %s\n", ss.str().c_str());

    vtkNew<vtkXMLStructuredGridWriter> writer;
    writer->SetFileName(ss.str().c_str());
    writer->SetCompressorTypeToZLib();
    writer->SetDataModeToBinary();
    writer->SetInputData(outData.GetPointer());
    writer->Write();

    return 0;
}
