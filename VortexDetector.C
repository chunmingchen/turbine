// An example detecting vortices in curvilinear grid data
// Samples regularly and compute the four criteria
// Stores into raw files with nrrd headers
// By Chun-Ming Chen

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <vtkXMLMultiBlockDataReader.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkInterpolatedVelocityField.h>
#include <vtkNew.h>
#include <vtkMultiBlockPLOT3DReader.h>
#include <system/path.h>
#include "MultiBlockVectorFieldVTK.h"
#include "OSUFlow.h"

vtkInterpolatedVelocityField *dum;

#define unit 0.005f

vtkSmartPointer<vtkMultiBlockDataSet> load_list(char *list_fname)
{
    //char list_filename[1024];
    //sprintf(list_filename, "%s/list", DATA_PATH);
    std::string path = getPath(list_fname);

    FILE *fp = fopen(list_fname, "rt");
    char s[1024];
    fgets(s, 1024, fp);

    int blocks = atoi(s);
    int i;
    vtkSmartPointer<vtkMultiBlockDataSet> mb = vtkMultiBlockDataSet::New();

    for (i=0; i<blocks; i++)
    {
        char file1[1024], file2[1024]; // q, xyz
        fgets(s, 1024, fp);
        *strchr(s, '\n')=0; // remove the last new-line
        sprintf(file1, "%s%s", path.c_str(), s);
        fgets(s, 1024, fp);
        *strchr(s, '\n')=0; // remove the last new-line
        sprintf(file2, "%s%s", path.c_str(), s);
        printf("xyz: [%s]   q: [%s]\n", file1, file2);

        // Start by loading some data.
        vtkNew<vtkMultiBlockPLOT3DReader> pl3dReader;
        pl3dReader->SetXYZFileName(file1);
        pl3dReader->SetQFileName(file2);
        pl3dReader->SetScalarFunctionNumber(100);
        pl3dReader->SetVectorFunctionNumber(200);
        pl3dReader->SetAutoDetectFormat(1);

#if 0
        pl3dReader->AddFunction(100);
        pl3dReader->AddFunction(110);
        pl3dReader->AddFunction(111);
        pl3dReader->AddFunction(112);
        pl3dReader->AddFunction(113);
        pl3dReader->AddFunction(120);
        pl3dReader->AddFunction(130);
        pl3dReader->AddFunction(140);
        pl3dReader->AddFunction(144);
        pl3dReader->AddFunction(153);
        pl3dReader->AddFunction(163);
        pl3dReader->AddFunction(170);
        pl3dReader->AddFunction(184);
        pl3dReader->AddFunction(211);
#endif

        pl3dReader->AddFunction(200);
        //pl3dReader->AddFunction(201);
        //pl3dReader->AddFunction(202);
        //pl3dReader->AddFunction(210);
        //pl3dReader->AddFunction(212);


        pl3dReader->Update();
        vtkDataSet *current_data = vtkDataSet::SafeDownCast(pl3dReader->GetOutput()->GetBlock(0));



        mb->SetBlock(i, current_data);

    }
    return mb;
}

int main(int argc, char ** argv)
{

	OSUFlow *osuflow = new OSUFlow();

	VECTOR3 minLen, maxLen, minB, maxB;
	minB[0] = 0; minB[1] = 0; minB[2] = 0;
    maxB[0] = 200; maxB[1] = 200; maxB[2] = 200;


    // load data
#if 0
    // vtk multiblock
	vtkXMLMultiBlockDataReader *reader = vtkXMLMultiBlockDataReader::New();
    reader->SetFileName(argv[1]);
    reader->Update();
    vtkSmartPointer<vtkMultiBlockDataSet> data = vtkMultiBlockDataSet::SafeDownCast(reader->GetOutput());
#else
    // plot 3d
    vtkSmartPointer<vtkMultiBlockDataSet> data = load_list(argv[1]);
#endif

	osuflow = new OSUFlow;
    CVectorField *field = new MultiBlockVectorFieldVTK( data.GetPointer() );
	osuflow->SetFlowField( field );
	osuflow->Boundary(minLen, maxLen);
	printf(" volume boundary X: [%f %f] Y: [%f %f] Z: [%f %f]\n",
								minLen[0], maxLen[0], minLen[1], maxLen[1],
								minLen[2], maxLen[2]);

	char out_fname[4][256];
	char prefix[4][10]={"lambda2", "q", "delta", "gamma"};
	int i;
	for (i=0; i<4; i++)
	{
		sprintf(out_fname[i], "%s_%s.raw", "vortex", prefix[i]);
		printf("Output file: %s\n", out_fname[i]);
	}

	int count=0;
	FILE *fp[4];
	for (i=0; i<4; i++)
		fp[i] = fopen(out_fname[i], "wb");

	float x,y,z;
	for (z=minLen[2]; z<=maxLen[2]; z+=unit) {
		for (y=minLen[1]; y<=maxLen[1]; y+=unit)
			for (x=minLen[0]; x<=maxLen[0]; x+=unit)
			{
#if 0
				VECTOR3 v;
				osuflow->GetFlowField()->at_phys(VECTOR3(x,y,z), 0, v);
				printf("%f %f %f\n", v[0], v[1], v[2]);
#endif
				float f[4]; //lambda2, q, delta, gamma;
				osuflow->GetFlowField()->GenerateVortexMetrics(VECTOR3(x,y,z), f[0], f[1], f[2], f[3], unit/2.f);
				for (i=0; i<4; i++)
					fwrite((char *)&f[i], 1, 4, fp[i]);
				count++;
			}
		printf("z=%f\n", z);
	}
	for (i=0; i<4; i++)
		fclose(fp[i]);


	// get dim for given range
	int bdim[3];
	bdim[0] = bdim[1] = bdim[2] = 0;
	{
		int d;
		for (d=0; d<3; d++)
			for (x=minLen[d]; x<=maxLen[d]; x+=unit)
				bdim[d] ++;
	}

	for (i=0; i<4; i++)
	{
		// get out_fname filename only (no path)
		char *out_fname_no_path = strrchr(out_fname[i], '/');
		if (out_fname_no_path==NULL) out_fname_no_path = out_fname[i]; else out_fname_no_path++;

		char out_desc_fname[256];
		sprintf(out_desc_fname, "%s_%s.nhdr", "vortex", prefix[i]);
		FILE *fp = fopen(out_desc_fname, "wt");
		fprintf(fp,
				"NRRD0004\n"
				"type: float\n"
				"dimension: 3\n"
				"sizes: %d %d %d\n"
				"encoding: raw\n"
				"data file: %s\n"
				"space origin: (%f,%f,%f)\n"
	            "space directions: (%f,0,0) (0,%f,0) (0,0,%f)\n"
				"# sampling distance: %f\n",
				bdim[0], bdim[1], bdim[2], out_fname_no_path, minLen[0], minLen[1], minLen[2], unit, unit , unit,
				unit);
		fclose(fp);
	}

	printf("Done (%d elems)\n", count);

	return 0;
}
