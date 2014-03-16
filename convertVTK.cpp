// Reference: VTK/Examples/VisualizationAlgorithms/Python/StreamlinesWithLineWidget.py
#include <stdio.h>
#include <stdlib.h>

#include <list>
#include <iterator>

//#include "vtkOSUFlow.h"
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
#include "vtkMultiBlockDataSet.h"
#include "vtkXMLMultiBlockDataWriter.h"
#include "vtkZLibDataCompressor.h"
// streamline
#include "vtkStreamLine.h"

using namespace std;

#define DATA_PATH "./" //"/data/flow2/Stg37/"

vtkLineWidget *lineWidget;
vtkRenderWindow *renWin;
vtkPolyData *seeds ;


vtkSmartPointer<vtkMultiBlockDataSet> load_list(char *list_fname)
{
	//char list_filename[1024];
	//sprintf(list_filename, "%s/list", DATA_PATH);
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
		sprintf(file1, "%s%s", DATA_PATH, s);
		fgets(s, 1024, fp);
		*strchr(s, '\n')=0; // remove the last new-line
		sprintf(file2, "%s%s", DATA_PATH, s);
		printf("xyz: [%s]   q: [%s]\n", file1, file2);

		// Start by loading some data.
		vtkNew<vtkMultiBlockPLOT3DReader> pl3dReader;
		pl3dReader->SetXYZFileName(file1);
		pl3dReader->SetQFileName(file2);
		pl3dReader->SetScalarFunctionNumber(100);
		pl3dReader->SetVectorFunctionNumber(200);
	    pl3dReader->SetAutoDetectFormat(1);

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

	    //pl3dReader->AddFunction(200);
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

int main(int argc, char **argv)
{
	printf("Usage: convertVTK <list file>\n");
	printf("Press 'i' to change the rake\n");

    // compressor
    vtkSmartPointer<vtkDataCompressor> compressor = vtkZLibDataCompressor::New();

	// write data
	vtkSmartPointer<vtkMultiBlockDataSet> mb = load_list(argv[1]);
	vtkNew<vtkXMLMultiBlockDataWriter> mbw;
	mbw->SetFileName("merged.vtm");
	mbw->SetDataModeToBinary();
	mbw->SetInputData(mb.GetPointer());
    mbw->SetCompressor(compressor);
	mbw->Write();


	return 0;
}



