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

struct Stat
{
	int nNegXVel;
	double in_massFlow, out_massFlow;
};

int countNegXVel(vtkSmartPointer<vtkImageData> data)
{
	vtkDataArray *array = data->GetPointData()->GetArray("Velocity");
	vtkDataArray *valid_array = data->GetPointData()->GetArray("vtkValidPointMask");
	int count = 0;
	int n = array->GetNumberOfTuples();
	for (int i=0; i<n; i++)
	{
		double valid;
		valid_array->GetTuple(i, &valid);
		if (valid != 0) {
			//printf("valid=%lf\n", valid);
			double value[3];
			array->GetTuple(i, value);
			if (value[0] < 0)
				count ++;
		}
		//rawDataStorage[i][z][y].push_back(value[vtiTupleIdx[i]]);
		//printf("%lg ", value[vtiTupleIdx[i]]);
	}
	return count;

}

void computeMassFlowRate(vtkSmartPointer<vtkImageData> data, double *in_rate, double *out_rate)
{
	vtkDataArray *varray = data->GetPointData()->GetArray("Velocity");
	vtkDataArray *darray = data->GetPointData()->GetArray("Density");
	vtkDataArray *valid_array = data->GetPointData()->GetArray("vtkValidPointMask");
	double sum_in = 0, sum_out = 0;

    int *ext = data->GetExtent();
    int dim[3];
    dim[0] = ext[1]+1;
    dim[1] = ext[3]+1;
    dim[2] = ext[5]+1;

	for (int z=0; z<dim[2]; z++)
		for (int y=0; y<dim[1]; y++)
			for (int in_out=0; in_out<2; in_out++)
			{
				int x;
				if (in_out==0) // in
					x=0;
				else
					x=dim[0]-1;

				int idx = x+dim[0]*(y+dim[1]*z);
				double valid;
				valid_array->GetTuple(idx, &valid);
				if (valid != 0) {
					double v[3];
					varray->GetTuple(idx, v);
					double vmag = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

					double d;
					darray->GetTuple(idx, &d);
					if (in_out==0)
						sum_in += d*vmag;
					else
						sum_out += d*vmag;
				}
			}
	*in_rate = sum_in;
	*out_rate = sum_out;
}

int main(int argc, char **argv)
{
	char *filename = argv[1];

	// load
    vtkSmartPointer<vtkXMLImageDataReader> reader = vtkXMLImageDataReader::New();
    reader->SetFileName(filename);
    reader->UpdateInformation();
    reader->Update();

    int *ext = reader->GetOutput()->GetExtent();
    fprintf(stderr, "Extent: %d %d %d %d %d %d\n", ext[0], ext[1], ext[2], ext[3], ext[4], ext[5]);

    // get data
    vtkSmartPointer<vtkImageData> data = reader->GetOutput() ;

    // get file id
    int id = atoi(argv[2]);

    Stat stat;
    stat.nNegXVel = countNegXVel(data);
    computeMassFlowRate(data, &stat.in_massFlow, &stat.out_massFlow );
    printf("%d, %d, %lg, %lg \n", id, stat.nNegXVel, stat.in_massFlow, stat.out_massFlow);
}
