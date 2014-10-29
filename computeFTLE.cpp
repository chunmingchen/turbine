#include <vector>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <OSUFlow.h>
#include <Field.h>


using namespace std;

int main(int argc, char **argv)
{
	printf("Usage: computeFTLE flowmap.raw w h d RES\n");
	int w = atoi(argv[2]),
		h = atoi(argv[3]),
		d = atoi(argv[4]);
	float RES = atof(argv[5]);
	printf("Opening file %s, w h d: %d %d %d\n", argv[1], w, h, d);

	vector<VECTOR3> offset(w*h*d);

	FILE *fp = fopen(argv[1], "rb");
	fread(&offset[0], w*h*d, 12, fp);
	fclose(fp);

	if (0) // now we generate flowmap
	{
		int x,y,z, count=0;
		for (z=0; z<d; z++)
			for (y=0; y<h; y++)
				for (x=0; x<w; x++) {
					offset[count] = offset[count] + VECTOR3(x*RES, y*RES, z*RES);
					//printf("offset: %f %f %f\n", offset[count][0], offset[count][1], offset[count][2]);
					count ++;
				}
	}

	// use osuflow
	OSUFlow *osuflow = new OSUFlow();
	float minB[3] = {0,0,0}, maxB[3];
	maxB[0] = w-1;  maxB[1] = h-1; maxB[2] = d-1;
	osuflow->CreateStaticFlowField((float *)&offset[0], w, h, d, minB, maxB);
	CVectorField *field = osuflow->GetFlowField();




	// FTLE
    vector<float> ftle(w*h*d);
	int x,y,z, count=0;
	for (z=0; z<d; z++)
	{
		for (y=0; y<h; y++)
			for (x=0; x<w; x++)
			{
                MATRIX3 jac = field->Jacobian(VECTOR3(x,y,z)),
                        jsquare; // TODO: delta in Jac. computation
				jac = jac * (1/RES);
				jsquare = jac.transpose() * jac;

				float m[3][3], eigenvalues[3];
				for (int i = 0; i < 3; i++) {
					for (int j = 0; j < 3; j++) {
						m[i][j] = jsquare(i, j);
					}
				}
				compute_eigenvalues(m, eigenvalues);
				float max_eig = max(eigenvalues[0], max(eigenvalues[1], eigenvalues[2]));
				//printf("%f\n", max_eig);
				ftle[count++] = log(max_eig)*.5;
			}
		printf("z=%d\n", z);
	}

	fp = fopen("ftle.raw", "wb");
	fwrite(&ftle[0], w*h*d, sizeof(float), fp);
	fclose(fp);

	printf("output: ftle.raw\n");
	return 0;
}
