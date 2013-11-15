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
#include "vtkXMLMultiBlockDataReader.h"
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkUnstructuredGrid.h"
// streamline
#include "vtkStreamLine.h"

using namespace std;

vtkLineWidget *lineWidget;
vtkOSUFlow *streamer;
vtkRenderWindow *renWin;
vtkPolyData *seeds ;

void computeStreamlines(vtkObject* caller, unsigned long eventId, void *clientdata, void *calldata)
{
	printf("compute\n");
	lineWidget->GetPolyData(seeds);
	renWin->Render();
	streamer->Update();
}

////////////////////










//////////////////

int main(int argc, char **argv)
{
	printf("Press 'i' to change the rake\n");

	// read data
	printf("Loading...\n");
#if 1
	vtkNew<vtkXMLMultiBlockDataReader> mbr;
	mbr->SetFileName("/data/flow2/Stg37/vtk/merged4.vtm");
	mbr->Update();
	vtkMultiBlockDataSet* mb = vtkMultiBlockDataSet::SafeDownCast( mbr->GetOutput() );

	vtkDataSet *data;
	data = vtkDataSet::SafeDownCast(mb->GetBlock(0));
#else
	vtkNew<vtkXMLUnstructuredGridReader> reader;
	reader->SetFileName("/data/flow2/Stg37/vtk/merged4.vtu");
	reader->Update();
	vtkUnstructuredGrid* data = vtkUnstructuredGrid::SafeDownCast( reader->GetOutput() );
#endif

	printf("Done loading\n");

	//
	// Determine seeds
	//
	// user can change the rake
	lineWidget = vtkLineWidget::New();
	lineWidget->SetInputData(data);
	lineWidget->SetResolution(21); // 22 seeds along the line
	lineWidget->SetAlignToYAxis();
	lineWidget->PlaceWidget();
	lineWidget->ClampToBoundsOn();
	seeds = vtkPolyData::New();
	lineWidget->GetPolyData(seeds);

#if 0
	//
	// vtkOSUFlow
	//
	streamer = vtkOSUFlow::New();
	streamer->SetInputData(data);
	streamer->SetSourceData(seeds);	//streamer->SetSourceConnection(rake->GetOutputPort());
	streamer->SetIntegrationDirectionToForward();
	streamer->SetMaximumPropagationTime(200);
	streamer->SetNumberOfThreads(1);
	streamer->VorticityOn();
	//streamer->getOSUFlow()->initOpenMP(8); // set number of processes to 8
#else
	// use vtkStreamLine
	vtkStreamLine *streamer = vtkStreamLine::New();
	streamer->SetInputData(data);
	//streamer->SetInputConnection(mbr->GetOutputPort());
	streamer->SetSourceData(seeds);
	streamer->SetMaximumPropagationTime(200);
	streamer->SetIntegrationStepLength(.001);
	streamer->SetStepLength(.001);
	streamer->SetNumberOfThreads(1);
	streamer->SetIntegrationDirectionToForward();//  IntegrateBothDirections();
	streamer->VorticityOff();
	printf("Streamer done\n");
#endif

	vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
	mapper->SetInputConnection(streamer->GetOutputPort());
	mapper->SetScalarRange(data->GetScalarRange());
	vtkActor *actor = vtkActor::New();
	actor->SetMapper(mapper);
	printf("Mapper done\n");


	//
	// outline
	//
	vtkStructuredGridOutlineFilter *outline = vtkStructuredGridOutlineFilter::New();
	outline->SetInputData(data);

	vtkPolyDataMapper *outlineMapper = vtkPolyDataMapper::New();
	outlineMapper->SetInputConnection(outline->GetOutputPort());

	vtkActor *outlineActor = vtkActor::New();
	outlineActor->SetMapper(outlineMapper);
	outlineActor->GetProperty()->SetColor(0,0,0);
	printf("Outline done\n");


	//
	// renderer
	//
	vtkRenderer *ren = vtkRenderer::New();
	renWin = vtkRenderWindow::New();
	renWin->AddRenderer(ren);
	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(renWin);
	vtkInteractorStyleTrackballCamera *style = vtkInteractorStyleTrackballCamera::New();
	iren->SetInteractorStyle(style);
	printf("renderer done\n");

	// line widget interactor
	lineWidget->SetInteractor(iren);
	lineWidget->SetDefaultRenderer(ren);
	vtkCallbackCommand *callback = vtkCallbackCommand::New();
	callback->SetCallback(computeStreamlines);
	lineWidget->AddObserver(vtkCommand::EndInteractionEvent, callback);
	printf("lineWidget done\n");

	//ren->AddActor(rakeActor);
	ren->AddActor(actor);
	ren->AddActor(outlineActor);
	ren->SetBackground(.5,.5,.5);

	renWin->SetSize(500,500);

	iren->Initialize();
	renWin->Render();
	iren->Start();

	return 0;
}



