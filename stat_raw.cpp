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

int dim[3];

#define xyz_id(x,y,z) ((x)+dim[0]*((y)+dim[1]*(z)))

// get tangent velocity at computation space x,y,z
VECTOR3 get_centroid_vec(vtkDataArray *ary, int offset[8][3] , int x, int y, int z)
{
    VECTOR3 vec[8];

    for (int off=0; off<8; off++)
    {
        // get velocity
        int idx = xyz_id(x+offset[off][0], y+offset[off][1], z+offset[off][2]);
        double *value = ary->GetTuple3(idx);
        vec[off] = VECTOR3(value[0], value[1], value[2]);
    }
    // centroid vec/pos
    VECTOR3 cvec = (vec[0]+vec[1]+vec[2]+vec[3]+vec[4]+vec[5]+vec[6]+vec[7])*0.125;
    return cvec;
}

VECTOR3 get_centroid_normal(vtkDataArray *ary, int offset[8][3], int x, int y, int z)
{
    VECTOR3 pos[8];

    for (int off=0; off<8; off++)
    {
        // get velocity
        int idx = xyz_id(x+offset[off][0], y+offset[off][1], z+offset[off][2]);
        double *value = ary->GetTuple3(idx);
        pos[off] = VECTOR3(value[0], value[1], value[2]);
    }
    VECTOR3 cpos0 = (pos[0] + pos[4])*0.5;
    VECTOR3 cpos1 = (pos[1] + pos[5])*0.5;
    VECTOR3 cpos2 = (pos[2] + pos[6])*0.5;
    VECTOR3 cpos3 = (pos[3] + pos[7])*0.5;

    // compute two vectors
    VECTOR3 va = cpos3-cpos0;
    VECTOR3 vb = cpos1-cpos2;

    VECTOR3 n = cross(vb, va);
    n.Normalize();
    return n;
}

#define DVT_VERBOSE
// Derivation of tangent velocity based on Prof. Chen
void computeDTangentVelocity(vtkSmartPointer<vtkMultiPieceDataSet> mb, char *out_filename)
{
    static int offset[8][3] = {{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0}, {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}};
    static int offset_2d[8][3] = {{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0}, {0,0,0}, {1,0,0}, {0,1,0}, {1,1,0}};
    vtkNew<vtkPoints> points;

    vtkNew<vtkCellArray> out_quad_array;

    vtkNew<vtkFloatArray> out_vel_array;
    out_vel_array->SetName("Velocity");
    out_vel_array->SetNumberOfComponents(3);

    vtkNew<vtkFloatArray> out_dvt_array;
    out_dvt_array->SetName("d_tangent_vel");
    out_dvt_array->SetNumberOfComponents(1);

    vtkNew<vtkFloatArray> out_normal_array;
    out_normal_array->SetName("normal");
    out_normal_array->SetNumberOfComponents(3);

#ifdef DVT_VERBOSE
    vtkNew<vtkFloatArray> out_vt1_array;
    out_vt1_array->SetName("Vt1");
    out_vt1_array->SetNumberOfComponents(1);
    vtkNew<vtkFloatArray> out_vt2_array;
    out_vt2_array->SetName("Vt2");
    out_vt2_array->SetNumberOfComponents(1);
    vtkNew<vtkFloatArray> out_vt3_array;
    out_vt3_array->SetName("Vt3");
    out_vt3_array->SetNumberOfComponents(1);
    vtkNew<vtkFloatArray> out_h2_array;
    out_h2_array->SetName("h2");
    out_h2_array->SetNumberOfComponents(1);
    vtkNew<vtkFloatArray> out_h3_array;
    out_h3_array->SetName("h3");
    out_h3_array->SetNumberOfComponents(1);
#endif

    // go through each piece:
    int pieces = mb->GetNumberOfPieces();
    for (int i=0; i<pieces; i++	)
    {
        //mb->GetPiece(i)->PrintSelf(std::cout, vtkIndent(0));
        vtkStructuredGrid *data = vtkStructuredGrid::SafeDownCast( mb->GetPiece(i) );
        assert(data);

        //< get extent
        int *ext = data->GetExtent();
        //int dim[3];
        dim[0] = ext[1]-ext[0]+1;
        dim[1] = ext[3]-ext[2]+1;
        dim[2] = ext[5]-ext[4]+1;
        //>

        vtkDataArray *v_ary = data->GetPointData()->GetArray("Velocity");
        vtkDataArray *p_ary = data->GetPoints()->GetData(); // point array
        int x,y,z; // data space

        std::vector<int> point_id_ary; //(dim[0]*dim[1]);

        {
            z = 0;
            for (y=0; y<dim[1]-1; y++)
                for (x=0; x<dim[0]-1; x++)
                {
                    VECTOR3 n1 = get_centroid_normal(p_ary, offset_2d, x,y,z);
                    VECTOR3 n2 = get_centroid_normal(p_ary, offset, x,y,z);
                    VECTOR3 n3 = get_centroid_normal(p_ary, offset, x,y,z+1);

                    VECTOR3 cvel1 = get_centroid_vec(v_ary, offset_2d, x,y,z);
                    VECTOR3 cvel2 = get_centroid_vec(v_ary, offset, x,y,z);
                    VECTOR3 cvel3 = get_centroid_vec(v_ary, offset, x,y,z+1);

                    // velocity in normal direction
                    float vn1 = dot(cvel1, n1);
                    float vn2 = dot(cvel2, n2);
                    float vn3 = dot(cvel3, n3);

                    // velocity in tangential direction
                    float vt1 = (cvel1-n1*vn1).GetMag();
                    float vt2 = (cvel2-n2*vn2).GetMag();
                    float vt3 = (cvel3-n3*vn3).GetMag();

                    VECTOR3 cpos1 = get_centroid_vec(p_ary, offset_2d, x,y,z);
                    VECTOR3 cpos2 = get_centroid_vec(p_ary, offset, x,y,z);
                    VECTOR3 cpos3 = get_centroid_vec(p_ary, offset, x,y,z+1);

                    float h2 = dot(cpos2-cpos1, n1);
                    float h3 = dot(cpos3-cpos1, n1);
                    double h2s = (double)h2*h2; //square
                    double h3s = (double)h3*h3;

                    float d_vt = ((h3s-h2s)*vt1 - h3s*vt2+h2s*vt3) / (h3*h2s-h2*h3s);

                    // VTK output :

                    // insert point
                    int id = points->InsertNextPoint(cpos2[0], cpos2[1], cpos2[2]);
                    point_id_ary.push_back(id);

                    // insert velocity
                    out_normal_array->InsertNextTuple3(n2[0], n2[1], n2[2]);
                    out_vel_array->InsertNextTuple3(cvel2[0], cvel2[1], cvel2[2]);
                    out_dvt_array->InsertNextTuple(&d_vt);

#ifdef DVT_VERBOSE
                    out_vt1_array->InsertNextTuple(&vt1);
                    out_vt2_array->InsertNextTuple(&vt2);
                    out_vt3_array->InsertNextTuple(&vt3);
                    out_h2_array->InsertNextTuple(&h2);
                    out_h3_array->InsertNextTuple(&h3);
#endif

                }
        }

        // create surface
        for (y=0; y<dim[1]-2; y++)
            for (x=0; x<dim[0]-2; x++)
            {

#define xyz_id1(x,y) ((x)+(dim[0]-1)*(y))
                vtkNew<vtkQuad> quad;
                quad->GetPointIds()->SetId(0, point_id_ary[xyz_id1(x,y)] );
                quad->GetPointIds()->SetId(1, point_id_ary[xyz_id1(x+1,y)] );
                quad->GetPointIds()->SetId(2, point_id_ary[xyz_id1(x+1,y+1)] );
                quad->GetPointIds()->SetId(3, point_id_ary[xyz_id1(x,y+1)] );
                out_quad_array->InsertNextCell(quad.GetPointer());
#undef xyz_id1
            }
    }
    printf("Done\n");
    vtkNew<vtkPolyData> poly;
    poly->SetPoints(points.GetPointer());
    poly->GetPointData()->SetVectors(out_vel_array.GetPointer());
    poly->GetPointData()->SetNormals(out_normal_array.GetPointer());
    poly->GetPointData()->SetScalars(out_dvt_array.GetPointer());
#ifdef DVT_VERBOSE
    poly->GetPointData()->AddArray(out_vt1_array.GetPointer());
    poly->GetPointData()->AddArray(out_vt2_array.GetPointer());
    poly->GetPointData()->AddArray(out_vt3_array.GetPointer());
    poly->GetPointData()->AddArray(out_h2_array.GetPointer());
    poly->GetPointData()->AddArray(out_h3_array.GetPointer());
#endif
    poly->SetPolys(out_quad_array.GetPointer());


    // save file
    printf("Output filename: %s\n", out_filename);
    vtkNew<vtkXMLPolyDataWriter> pw;
    pw->SetFileName(out_filename);
    pw->SetInputData(poly.GetPointer());
    pw->Write();

}


void computeParVel_simple(vtkSmartPointer<vtkMultiPieceDataSet> mb, int layer_id, char *out_filename)
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
        //int dim[3];
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
					double projected_vel =
							v1.GetMag()? dot(v1, vel)/v1.GetMag() : 0;
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
                quad->GetPointIds()->SetId(0, point_id_ary[xyz_id(x,y,0)] );
                quad->GetPointIds()->SetId(1, point_id_ary[xyz_id(x+1,y,0)] );
                quad->GetPointIds()->SetId(2, point_id_ary[xyz_id(x+1,y+1,0)] );
                quad->GetPointIds()->SetId(3, point_id_ary[xyz_id(x,y+1,0)] );
				out_quad_array->InsertNextCell(quad.GetPointer());
			}
	}
	vtkNew<vtkPolyData> poly;
	poly->SetPoints(points.GetPointer());
	poly->GetPointData()->SetVectors(out_vel_array.GetPointer());
	poly->GetPointData()->SetScalars(out_surface_vel_array.GetPointer());
	poly->SetPolys(out_quad_array.GetPointer());


	// save file
	vtkNew<vtkXMLPolyDataWriter> pw;
	pw->SetFileName(out_filename);
	pw->SetInputData(poly.GetPointer());
	pw->Write();

}

int main(int argc, char **argv)
{
    printf("usage: input_filename  output_filename\n");
    char *filename = argv[1];
    char *out_filename = argv[2];

	// load data
	vtkSmartPointer<vtkMultiPieceDataSet> mb = load_list(filename);


	printf("computing...\n");
    computeDTangentVelocity(mb, out_filename);


}
