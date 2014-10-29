
#include <vtkVersion.h>
#include <vtkXMLMultiBlockDataReader.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkCompositeDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkDataSet.h>
#include <vtkDataObject.h>
#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkFieldData.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkInformation.h>
#include <vtkContourFilter.h>
#include <vtkCompositeDataGeometryFilter.h>
#include <vtkProperty.h>
#include <vtkStructuredGridOutlineFilter.h>
#include <vtkTransform.h>
#include <vtkAxesActor.h>
#include <vtkCompositePolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkSliderRepresentation2D.h>
#include <vtkSliderWidget.h>
#include <vtkProperty2D.h>
#include <vtkCommand.h>
#include <vtkCallbackCommand.h>
#include <vtkTextProperty.h>
#include <vtkNew.h>
#include <vtkMultiBlockPLOT3DReader.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkZLibDataCompressor.h>
#include <vtkSphereSource.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkDoubleArray.h>
#include <vtkCamera.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkOutlineFilter.h>
#include <vtkLookupTable.h>
#include <vtkFloatArray.h>
#include <vtkVersion.h>
#include <vtkPoints.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataReader.h>
#include <vtkMultiPieceDataSet.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkQuad.h>
#include <vtkErrorCode.h>
#include <vtkThreshold.h> // for vortex angle statistics
#include <vtkConnectivityFilter.h> // for vortex angle statistics
#include <vtkXMLUnstructuredGridWriter.h> // for testing
#include <vtkUnstructuredGrid.h>
#include <vtkMath.h>
#include <vtkExtractGrid.h>

//Definition of methods
int compute_3d_to_1d_map(int x,int y,int z, int dimx, int dimy, int dimz)
{
    return x + dimx*(y+dimy*z);
}

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

    for (i=0; i<blocks; i++)
    {
        char file1[1024], file2[1024]; // q, xyz
        fgets(s, 1024, fp);
        *strchr(s, '\n')=0; // remove the last new-line
        sprintf(file1, "%s/%s", DATA_PATH, s);
        fgets(s, 1024, fp);
        *strchr(s, '\n')=0; // remove the last new-line
        sprintf(file2, "%s/%s", DATA_PATH, s);
        printf("xyz: [%s]   q: [%s]\n", file1, file2);

        // Start by loading some data.
        vtkNew<vtkMultiBlockPLOT3DReader> reader;
        reader->SetXYZFileName(file1);
        reader->SetQFileName(file2);
        reader->SetScalarFunctionNumber(110);
        reader->SetVectorFunctionNumber(202);
        reader->SetAutoDetectFormat(1);

        //reader->AddFunction(100); //density
        reader->AddFunction(110); //pressure
        //reader->AddFunction(120); //temp
        //reader->AddFunction(130); //enthalpy
        //reader->AddFunction(140); //internal energy
        //reader->AddFunction(144); //kinetic energy
        //reader->AddFunction(153); //vel magnitude
        //reader->AddFunction(163); //stagnation energy
        //reader->AddFunction(170); //entropy
        //reader->AddFunction(184); //swirl
        //reader->AddFunction(211); //vorticity magnitude

        //available vector fields in the data
        //reader->AddFunction(200); //velocity
        //reader->AddFunction(201); //vorticity
        reader->AddFunction(202); //momentum
        //reader->AddFunction(210); //pressure gradient
        //reader->AddFunction(212); //starin rate

        reader->Update();
        vtkDataSet *current_data = vtkDataSet::SafeDownCast(reader->GetOutput()->GetBlock(0));

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

        mb->SetPiece(i, current_data);


    }
    return mb;
}


int main()
{
    //int numBlocks = mb->GetNumberOfBlocks();
    double current_block_bounds[6];
    int dim[3] = {151,71,56};
    int idx=0;
    string fpath = DATA_DIR;
    double pt[3];
    int startT = 0;
    int endT = 50;
    int trange = endT - startT + 1;
    vector<vector<float> > pointAry(trange, vector<float>(3, 0));
    vector<vector<float> > pressureAry;

    //vtkSmartPointer<vtkCompositeDataProbeFilter> prober = vtkCompositeDataProbeFilter::New();

    //for(int bb=21;bb<=30;bb++)
    //{
    //Iterate over all the blocks to extract the physical locations for probing
    //vtkSmartPointer<vtkPoints> probePoints = vtkSmartPointer<vtkPoints>::New();
    //vtkDataSet *current_block = vtkDataSet::SafeDownCast(mb->GetBlock(0));
    //current_block->GetBounds(current_block_bounds);

    int i=0;
    for(int z=0;z<56;z+=14) //dim[2] = z
    {
        for(int y=40;y<60;y+=5) //dim[1] = y
        {
            for(int x=40;x<60;x+=5) //dim[0] = x
            {
                //idx = compute_3d_to_1d_map(x,y,z,dim[0],dim[1],dim[2]);
                //current_block->GetPoint(idx,pt);
                //probePoints->InsertNextPoint(pt);
                pointAry[i][0] = x;
                pointAry[i][1] = y;
                pointAry[i][2] = z;
                i++;
            }
        }
    }

    int t;
    for (t=startT; t<=endT; t++)
    {
        int b;
        char fname[1024];
        sprintf(fname, "%s/%d.list", DATA_DIR, t+1);
        vtkSmartPointer<vtkMultiPieceDataSet> mb = load_list(fname);
        for (b=0; b<36; b++)
        {
            vtkDataSet *current_block = vtkDataSet::SafeDownCast(mb->GetBlock(b));

        }
        //current_block->GetBounds(current_block_bounds);

    }

    //vtkSmartPointer<vtkPolyData> probePolyData = vtkSmartPointer<vtkPolyData>::New();
    //probePolyData->SetPoints(probePoints);

    int totPtsEachProbe = 4*4*4;
    float **probedPresVals;

    probedPresVals = (float **)malloc(totPtsEachProbe*sizeof(float *));
    for(int jj=0;jj<totPtsEachProbe;jj++)
        probedPresVals[jj] = (float *)malloc(trange*sizeof(float));

    for(int ii=0;ii<totPtsEachProbe;ii++)
        for(int jj=0;jj<trange;jj++)
            probedPresVals[ii][jj] = 0.0;
    //create outputfile
    string fileName;
    stringstream blk;
    blk<<probedBlock;
    fileName = fpath + "PressureProbeBlock" + blk.str() + ".txt";
    ofstream fp;
    fp.open(fileName.c_str());
    int blockToLoad=36; //change it to 36 TODO
    int scalarVarList[14] = {1,1,0,0,0,0,0,0,0,0,0,0,0,0};

    //Load each time step and probe pressure field to get the values
    for(int i=startT;i<=endT;i++)
    {
        //load each time step
        int tt = 6201 + (i-1)*25;
        generate_list(tt,blockToLoad);
        string fname = "input_list.txt";
        vtkSmartPointer<vtkMultiBlockDataSet> mbdata = load_data_listFile(fname,scalarVarList);

        //Create probe filter and extract the pressure values
        prober->SetSourceData(mbdata);
        prober->SetInputData(probePolyData);
        prober->SpatialMatchOn();
        prober->Update();
        //prober->GetOutput()->GetPointData()->PrintSelf(cout,vtkIndent(0));

        //update the probed values
        vtkFloatArray* pressureData = vtkFloatArray::SafeDownCast (prober->GetOutput()->GetPointData()->GetArray("Pressure"));

        /*vtkSmartPointer<vtkPolyData> polydata = prober->GetPolyDataOutput();
        cout<<polydata->GetNumberOfPoints()<<endl;

        vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        writer->SetInputData(polydata);
        writer->SetFileName("x.vtp");
        writer->Write();*/

        for(int k=0;k<pressureData->GetNumberOfTuples();k++)
            probedPresVals[k][i-startT] = pressureData->GetValue(k);
    }

    fp<<trange<<endl;

    //Write the output file for the prober
    for(int i=0;i<totPtsEachProbe;i++)
    {
        double *point = probePoints->GetPoint(i);

        fp<<point[0]<<" "<<point[1]<<" "<<point[2]<<" ";

        for(int j=startT;j<=endT;j++)
        {
            fp<< probedPresVals[i][j-startT]<<" ";
        }

        fp<<endl;
    }

    fp.close();

}
