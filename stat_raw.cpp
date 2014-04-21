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
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkQuad.h>

#include "VectorMatrix.h"


vtkSmartPointer<vtkMultiPieceDataSet> load_list(char *list_fname)
{
	//char list_filename[1024];
	//sprintf(list_filename, "%s/list", DATA_PATH);
	FILE *fp = fopen(list_fname, "rt");
	char s[1024];
	fgets(s, 1024, fp);

	int blocks = atoi(s);
	int i;
	vtkSmartPointer<vtkMultiPieceDataSet> mb = vtkMultiPieceDataSet::New();

	std::string path = getPath(list_fname);

	for (i=0; i<blocks; i++)
	{
		char file1[1024], file2[1024]; // q, xyz
		fgets(s, 1024, fp);
		*strchr(s, '\n')=0; // remove the last new-line
		sprintf(file1, "%s/%s", path.c_str(), s);
		fgets(s, 1024, fp);
		*strchr(s, '\n')=0; // remove the last new-line
		sprintf(file2, "%s/%s", path.c_str(), s);
		printf("xyz: [%s]   q: [%s]\n", file1, file2);

		// Start by loading some data.
		vtkNew<vtkMultiBlockPLOT3DReader> reader;
		reader->SetXYZFileName(file1);
		reader->SetQFileName(file2);
		reader->SetScalarFunctionNumber(110);
		reader->SetVectorFunctionNumber(200);
		reader->SetAutoDetectFormat(1);

	    //reader->AddFunction(100); //density
	    reader->AddFunction(110); //pressure
	    //reader->AddFunction(120); //temp
	    ////reader->AddFunction(130); //enthalpy
	    ////reader->AddFunction(140); //internal energy
	    ////reader->AddFunction(144); //kinetic energy
	    ////reader->AddFunction(153); //vel magnitude
	    ////reader->AddFunction(163); //stagnation energy
	    //reader->AddFunction(170); //entropy
	    ////reader->AddFunction(184); //swirl
	    //reader->AddFunction(211); //vorticity magnitude

	    //available vector fields in the data
	    reader->AddFunction(200); //velocity
	    ////reader->AddFunction(201); //vorticity
	    ////reader->AddFunction(202); //momentum
	    ////reader->AddFunction(210); //pressure gradient
	    ////reader->AddFunction(212); //starin rate

	    reader->Update();
		vtkDataSet *current_data = vtkDataSet::SafeDownCast(reader->GetOutput()->GetBlock(0));

#if 0
		//extract uvel
		//<
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
	    //>
#endif

	    mb->SetPiece(i, current_data);


	}
	return mb;
}


void computeParVel(vtkSmartPointer<vtkMultiPieceDataSet> mb, int layer_id, char *out_filename)
{
	vtkNew<vtkPoints> points;

	// go through each piece:
	int pieces = mb->GetNumberOfPieces();
	int i;
	//pieces=1;

	vtkNew<vtkFloatArray> out_vel_array;
	out_vel_array->SetName("Velocity");
	out_vel_array->SetNumberOfComponents(3);

	vtkNew<vtkCellArray> out_quad_array;

	vtkNew<vtkFloatArray> out_surface_vel_array;
	out_surface_vel_array->SetName("Surface_vel");
	out_surface_vel_array->SetNumberOfComponents(1);

	for (i=0; i<pieces; i++	)
	{
		//mb->GetPiece(i)->PrintSelf(std::cout, vtkIndent(0));
		vtkStructuredGrid *data = vtkStructuredGrid::SafeDownCast( mb->GetPiece(i) );
		assert(data);

		//< get extent
	    int *ext = data->GetExtent();
	    int dim[3];
	    dim[0] = ext[1]+1;
	    dim[1] = ext[3]+1;
	    dim[2] = ext[5]+1;
	    assert(ext[0]==0);
	    assert(ext[2]==0);
	    assert(ext[4]==0);
	    //>

		vtkDataArray *vel_array = data->GetPointData()->GetArray("Velocity");
		int x,y,z; // data space
		z = layer_id;

		std::vector<int> point_id_ary(dim[0]*dim[1]);

#define xyz_id(x,y,z) (x+dim[0]*(y+dim[1]*z))
		for (y=0; y<dim[1]; y++)
			for (x=0; x<dim[0]; x++)
			{
				int idx = xyz_id(x,y,z);

				// get physical space
				double p[3]; // physical space
				data->GetPoint(x,y,z, p);

				// insert point
				int id = points->InsertNextPoint(p);
				point_id_ary.push_back(id);

				// get velocity
				double *value = vel_array->GetTuple3(idx);
				out_vel_array->InsertNextTuple(value);

				// compute out surface velocity
				if (x<dim[0]-1) {
					double p1[3];
					data->GetPoint(x+1,y,z, p1);

					VECTOR3 v1(p1[0]-p[0], p1[1]-p[1], p1[2]-p[2]);
					VECTOR3 vel(value[0], value[1], value[2]);
					double projected_vel = dot(v1, vel)/v1.GetMag();
					out_surface_vel_array->InsertNextTuple(&projected_vel);
				} else {
					double v = 0;
					out_surface_vel_array->InsertNextTuple(&v);
				}




			}

		// create surface
		for (y=0; y<dim[1]-1; y++)
			for (x=0; x<dim[0]-1; x++)
			{
				vtkNew<vtkQuad> quad;
				quad->GetPointIds()->SetId(0, point_id_ary[xyz_id(x,y,z)] );
				quad->GetPointIds()->SetId(1, point_id_ary[xyz_id(x+1,y,z)] );
				quad->GetPointIds()->SetId(2, point_id_ary[xyz_id(x+1,y+1,z)] );
				quad->GetPointIds()->SetId(3, point_id_ary[xyz_id(x,y+1,z)] );
				out_quad_array->InsertNextCell(quad.GetPointer());
			}
	}
	vtkNew<vtkPolyData> poly;
	poly->SetPoints(points.GetPointer());
	poly->GetPointData()->SetVectors(out_vel_array.GetPointer());
	poly->GetPointData()->SetScalars(out_surface_vel_array.GetPointer());
	poly->SetPolys(out_quad_array.GetPointer());


	vtkNew<vtkXMLPolyDataWriter> pw;
	pw->SetFileName(out_filename);
	pw->SetInputData(poly.GetPointer());
	pw->Write();

}

int main(int argc, char **argv)
{
	printf("usage: filename layer_id output_filename\n");
	char *filename = argv[1];
	int layer_id = atoi(argv[2]);
	char *out_filename = argv[3];

	// load data
	vtkSmartPointer<vtkMultiPieceDataSet> mb = load_list(filename);


	printf("computing...\n");
    computeParVel(mb, layer_id, out_filename);


}
