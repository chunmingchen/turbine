// Reference: VTK/Examples/VisualizationAlgorithms/Python/StreamlinesWithLineWidget.py
#include <stdio.h>
#include <stdlib.h>

#include <list>
#include <iterator>

#include "vtkOSUFlow.h"
#include "vtkDataSet.h"
#include "vtkStructuredGrid.h"
#include "vtkMultiBlockPLOT3DReader.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkSmartPointer.h"
#include "vtkLineSource.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkStructuredGridOutlineFilter.h"
#include "vtkProperty.h"
#include "vtkLineWidget.h"
#include "vtkCommand.h"
#include "vtkCallbackCommand.h"
#include "vtkNew.h"
#include "vtkSmartPointer.h"
#include "vtkStructuredGrid.h"
#include "vtkUnstructuredGrid.h"
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkPointData.h"
#include "vtkCellLocatorInterpolatedVelocityField.h"
// streamline
#include "vtkStreamLine.h"
#include "VectorFieldVTK.h"
#include "OSUFlow.h"

using namespace std;

#define DATA_PATH "/data/flow2/Stg37/"

vtkLineWidget *lineWidget;
vtkOSUFlow *streamer;
vtkRenderWindow *renWin;
vtkPolyData *seeds ;

vtkSmartPointer<vtkUnstructuredGrid> getData_vtu()
{
	vtkNew<vtkXMLUnstructuredGridReader> reader;
	reader->SetFileName("/data/flow2/Stg37/vtk/merged4.vtu");
	reader->Update();

	vtkUnstructuredGrid *data = reader->GetOutput();
	if (data==0) {
		printf("data error\n");
	}
	data->Print(cout);

	data->GetPointData()->Print(cout);
	int r = data->GetPointData()->SetActiveAttribute("Velocity", vtkDataSetAttributes::VECTORS);
	printf("index=%d\n", r);

#if 0
	//reader->Delete();
	// test
	vtkNew<vtkCellLocatorInterpolatedVelocityField> interpolator;
	interpolator->AddDataSet(data);
	double vel[3];
	double coords[4] = {0,0,0,0};
	for (int k=0; k<3; k++)
	{
		bool success = (interpolator)->FunctionValues(coords, vel);
		printf("querry success = %d, vec = %lf %lf %lf\n", success, vel[0], vel[1], vel[2]);
	}
#endif

	return data;
}

#define RES 0.002 // sampling resolution
int main(int argc, char **argv)
{
	// read data
	printf("Loading...\n");
	vtkSmartPointer<vtkUnstructuredGrid> data = getData_vtu();
	double *bounds = data->GetBounds();
	printf("Boundary: %lf %lf %lf %lf %lf %lf\n", bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);

	int w,h,d;
	w = (bounds[1]-bounds[0])/RES+1;
	h = (bounds[3]-bounds[2])/RES+1;
	d = (bounds[5]-bounds[4])/RES+1;
	printf("w=%d h=%d d=%d\n", w,h,d);


	int concurrent_seeds = w*h;
	vector<VECTOR3> seedAry(concurrent_seeds);

	OSUFlow *osuflow = new OSUFlow();
	osuflow->SetFlowField(new VectorFieldVTK(data));
	osuflow->SetIntegrationParams(.001, .01);

	FILE *fp = fopen("flowmap.raw", "wb");
	vector<VECTOR3> endPos(concurrent_seeds);

	int x,y,z, i;
//#pragma omp for
	for (z=0; z<d; z++)
	{
		seedAry.clear();
		i=0;
		for (y=0; y<h; y++)
			for (x=0; x<w; x++)
				seedAry[i++] = VECTOR3(x*RES+bounds[0],y*RES+bounds[2],z*RES+bounds[4]);

		int cur_seeds = concurrent_seeds;
		//if (cur_seeds+i > seeds) cur_seeds = seeds-i;

		i=0;
		// gen
		list<vtListSeedTrace*> trace;
		osuflow->GenStreamLines(&seedAry[i], FORWARD_DIR, cur_seeds, 1000, trace);
		if (trace.size() != cur_seeds)
			printf("Error: list gets %zu seeds != %d", trace.size(), cur_seeds);


		int j = 0;
		std::list<vtListSeedTrace*>::iterator it;
		for (it = trace.begin(); it != trace.end(); it++)
		{
			vtListSeedTrace::iterator it1 ;
			//printf("trace: %zu \n ", (*it)->size());
			if ((*it)->size()) {

				it1 = (*it)->end();
				--it1;

				VECTOR3 &p = **it1;
				//VECTOR3 &p0 = seedAry[i+j];
				endPos[j] = p; //VECTOR3(p[0]-p0[0], p[1]-p0[1], p[2]-p0[2]);
				//printf("%f %f %f\n", p[0], p[1], p[2]);

			} else {
				VECTOR3 &p0 = seedAry[i+j];
				endPos[j] = p0;
			}

			for (it1=(*it)->begin(); it1!=(*it)->end(); it1++)
				delete (*it1);
			delete (*it);
			j++;
		}

		fwrite(&endPos[0], 12, cur_seeds, fp);
		printf("z=%d done\n", z);
	}


	delete osuflow;

	fclose(fp);
	printf("w=%d h=%d d=%d RES=%f\n", w,h,d, RES);

		// init results
	printf(" done integration.  File: flowmap.raw\n");



	return 0;
}



