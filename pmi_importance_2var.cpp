///////////////////////////////////////
// Soumya Dutta                      //  
// OSU CSE                           //
// Spring 2014                       //
///////////////////////////////////////
//Compile: g++ pmi_importance_2var.cpp -o pmi_importance_2var -I/home/soumya/
//Run: ./pmi_importance_2var 250 250 50 2 256 25

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

//GLM includes
#include <glm/glm.hpp>

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



#define pi 3.14

using namespace std;
using namespace glm;

int xdim;
int ydim;
int zdim;
int selectData;
int N;
int timestep=0;

int ARR_DIM;
int *Array1; 
int *Array2; 
int **ArrayComb; 
float *maxValsFirst;
float *minValsFirst;
float *maxValsSecond;
float *minValsSecond;
vector < vector <float> > PMI_Vals;
vector < vector <float> > Joint_Dist_Vals;
float **CosineSimilariy;
float dataMaxFirst,dataMinFirst,dataMaxSecond,dataMinSecond;
vector < vector <vector <vector <float> > > > rawDataStorage;
const int numFiles = 2;
string variable_name;
float *temp_value;
float delta=0.00000032;

float *I11,*I12,*I21,*I22,*I31,*I32;

map<float,vec2> pairMap;

// VTI
bool isVTI = 0;
char *vtiVarname[2];
int vtiTupleIdx[2]={0,1}; // get which element of the array ?

float log2(float x)
{
	return log(x)/log(2.0);
}
//This function computes all the I-metrices
void compute_metrics()
{
	float prob_of_x_given_y=0.0, prob_of_y_given_x=0.0, prob_of_x=0.0, prob_of_y=0.0;
	int x=0;

	for(int i=0;i<ARR_DIM;i++)
	{
		for(int j=0;j<ARR_DIM;j++)
		{

			if(Array1[i]==0)
			{
				prob_of_y_given_x = 0;
			}
			else
			{
				prob_of_y_given_x = (float)ArrayComb[i][j] / (float)Array1[i];
			}	

			prob_of_y = (float)Array2[j] / (float)N;

			if(prob_of_y_given_x != 0 && prob_of_y != 0)
			{
				I11[i] += prob_of_y_given_x * log(prob_of_y_given_x / prob_of_y);
			}
			if(prob_of_y_given_x != 0)
			{
				I21[i] += prob_of_y_given_x * log(prob_of_y_given_x);// - prob_of_y * log(prob_of_y);
			}
			if(prob_of_y != 0)
			{
				I21[i] -=  prob_of_y * log(prob_of_y);
			}

			if(Array2[i] == 0)
			{
				prob_of_x_given_y = 0;
			}
			else
			{
				//	prob_of_x_given_y = (float)ArrayComb[i][j] / (float)Array2[i]; // Array2[j]
				prob_of_x_given_y = (float)ArrayComb[j][i] / (float)Array2[i]; // Array2[j]
			}

			prob_of_x = (float)Array1[j] / (float)N;
			if(prob_of_x_given_y != 0 && prob_of_x != 0)
			{
				I12[i] += prob_of_x_given_y * log(prob_of_x_given_y / prob_of_x);
			}
			if(prob_of_x_given_y != 0)
			{
				I22[i] += prob_of_x_given_y * log(prob_of_x_given_y);// - prob_of_x * log(prob_of_x);
			}
			if(prob_of_x != 0)
			{
				I22[i] -= prob_of_x * log(prob_of_x);
			}

			if(prob_of_y_given_x > 1.0)
			{
				cout<<"Ooopps. value of prob_of_x_given_y is "<<prob_of_y_given_x<<" for i="<<i<<" and j="<<j<<endl;
				getchar();
			}

			if(prob_of_x_given_y > 1.0)
			{
				cout<<"Ooopps. value of prob_of_x_given_y is "<<prob_of_x_given_y<<" for i="<<i<<" and j="<<j<<endl;
				getchar();
			}


		}

	}

	//New Added
	for(int i=0;i<ARR_DIM;i++)
	{
		for(int j=0;j<ARR_DIM;j++)
		{
			if(Array1[i] == 0)
			{
				prob_of_y_given_x = 0;
			}
			else
			{
				prob_of_y_given_x = (float)ArrayComb[i][j] / (float)Array1[i];
			}

			prob_of_y = (float)Array2[j] / (float)N;
			I31[i] += prob_of_y_given_x * I22[j];

			if(Array2[i] == 0)
			{
				prob_of_x_given_y = 0;
			}
			else
			{
				prob_of_x_given_y = (float)ArrayComb[j][i] / (float)Array2[i]; // 
			}
			prob_of_x = (float)Array1[j] / (float)N;
			I32[i] += prob_of_x_given_y * I21[j];		
		}
	}
}

//construct a filename
string constructFileName(int varid, int j)
{
  string var_name;
  string fileName;
  string path;    
  string lastPart;
  string zeroes;

  //Combustion Data
  if(selectData==0)
  {
    if(varid==1)
      var_name = "jet_chi_";
    else if(varid==2)
      var_name = "jet_mixfrac_";
    else if(varid==3)
      var_name = "jet_hr_";
    else if(varid==4)
      var_name = "jet_vort_";
    else if(varid==5)
      var_name = "jet_Y_OH_";

    lastPart = ".dat_2_subsampled";
    path = "/home/soumya/Test_DataSet/Combustion/All_raw_subsampled_2/";

    //Construct the filename as per the value of current j
    stringstream tempstream2;
    if(j<10)
    {
      zeroes = "000";
      path = path + var_name + zeroes;
      tempstream2 << j;

    }
    else if(j>=10 && j<100)
    {
      zeroes = "00";
      path = path + var_name + zeroes;
      tempstream2 << j;
    }
    else
    {
      zeroes = "0";
      path = path + var_name + zeroes;
      tempstream2 << j;
    }

    //Construct the second file name
    fileName = path + tempstream2.str() + lastPart;
  }
  //Ion Front Data
  else if (selectData==1)
  {

    string ion_name = "multifield.";    

    if(varid==1)
      var_name = ".txt.density.f";
    else if(varid==2)
      var_name = ".txt.temp";         
    else if(varid==3)
      var_name = ".txt.H";
    else if(varid==4)
      var_name = ".txt.H2P";
    else if(varid==5)
      var_name = ".txt.H2";
    else if(varid==6)
      var_name = ".txt.He";         
    else if(varid==7)
      var_name = ".txt.HeP";
    else if(varid==8)
      var_name = ".txt.HePP";
    else if(varid==9)
      var_name = ".txt.HM"; 
    else if(varid==10)
      var_name = ".txt.HP"; 
    else if(varid==11)
      var_name = ".txt.curlmag"; 


    //lastPart = ".raw_normalized_2_subsampled";
    //path = "/media/jimmy/soumya/Test_DataSet/Ion_Front/All_timesteps_2/";

    lastPart = ".raw_2_subsampled"; 
    path = "/home/soumya/Test_DataSet/Ion_Front/All_raw_timesteps_2/"; 


    stringstream tempstream2;
    if(j<10)
    {
      zeroes = "000";
      tempstream2 << j;
      path = path + ion_name + zeroes + tempstream2.str() + var_name ;


    }
    else if(j>=10 && j<100)
    {
      zeroes = "00";
      tempstream2 << j;
      path = path + ion_name + zeroes + tempstream2.str() + var_name ;

    }
    else
    {
      zeroes = "0";
      tempstream2 << j;
      path = path + ion_name + zeroes + tempstream2.str() + var_name;

    }

    //Construct the second file name
    fileName = path + lastPart;
  }
  //Isabel Data
  else if (selectData==2)
  {
    if(varid==1)
      var_name = "Pf";
    else if(varid==2)
      var_name = "CLOUDf";          
    else if(varid==3)
      var_name = "PRECIPf";
    else if(varid==4)
      var_name = "QCLOUDf";
    else if(varid==5)
      var_name = "QGRAUPf";
    else if(varid==6)
      var_name = "QICEf";         
    else if(varid==7)
      var_name = "QRAINf";
    else if(varid==8)
      var_name = "QSNOWf";
    else if(varid==9)
      var_name = "QVAPORf"; 
    else if(varid==10)
      var_name = "TCf";
    else if(varid==11)
      var_name = "Uf";
    else if(varid==12)
      var_name = "Vf";
    else if(varid==13)
      var_name = "Wf";

    //lastPart = ".binLE.raw_corrected_normalized_2_subsampled";
    //path = "/media/jimmy/soumya/Test_DataSet/Isabel/All_normalized_timesteps_2/";  

    lastPart = ".binLE.raw_corrected_2_subsampled"; 
    path = "/home/soumya/Test_DataSet/Isabel/All_raw_subsampled_2/";


    //Construct the filename as per the value of current j
    stringstream tempstream2;
    if(j<10)
    {
      zeroes = "0";
      path = path + var_name + zeroes;
      tempstream2 << j;
    }
    else if(j>=10)
    {
      path = path + var_name;
      tempstream2 << j;
    }   

    //Construct the second file name
    fileName = path + tempstream2.str() + lastPart;   
  }
  //POP Data
  else if (selectData==3)
  {
    if(varid==1)
      var_name = "UVEL_";
    else if(varid==2)
      var_name = "VVEL_";         
    else if(varid==3)
      var_name = "TEMP_";
    else if(varid==4)
      var_name = "SALT_";
    else if(varid==5)
      var_name = "KVMIX_M_";
    else if(varid==6)
      var_name = "VDC_T_";          
    else if(varid==7)
      var_name = "VDC_S_";
    else if(varid==8)
      var_name = "RHO_";
    else if(varid==9)
      var_name = "UET_";  
    else if(varid==10)
      var_name = "VNT_";
    else if(varid==11)
      var_name = "WTT_";
    else if(varid==12)
      var_name = "UES_";
    else if(varid==13)
      var_name = "VNS_";
    else if(varid==14)
      var_name = "WTS_";
    else if(varid==15)
      var_name = "Q_";


    lastPart = ".raw_corrected_normalized_2_subsampled";
    path = "/home/soumya/Test_DataSet/LANL_POP_Yearly_Avg/All_timesteps_2/";

    //Construct the filename as per the value of current j
    stringstream tempstream2;
    path = path + var_name; 
    tempstream2 << j; 
    fileName = path + tempstream2.str() + lastPart;
  }

  if(selectData==4)
  {
  	if(varid==1)
      var_name = "Velocity_combined_time";
    else if(varid==2)
      var_name = "Pressure_combined_time";
    else if(varid==3)
      var_name = "Density_combined_time";
    else if(varid==4)
      var_name = "Temperature_combined_time";
    else if(varid==5)
      var_name = "Enthalpy_combined_time";
    else if(varid==6)
      var_name = "KineticEnergy_combined_time";
    else if(varid==7)
      var_name = "VelocityMagnitude_combined_time";
    else if(varid==8)
      var_name = "StagnationEnergy_combined_time";
    else if(varid==9)
      var_name = "Entropy_combined_time";
    else if(varid==10)
      var_name = "Swirl_combined_time";
    else if(varid==11)
      var_name = "VorticityMagnitude_combined_time";

  	lastPart = ".raw";
    path = "/home/soumya/Test_DataSet/Turbine/All_raw/";

    stringstream tempstream2;
    path = path + var_name; 
    tempstream2 << j; 
    fileName = path + tempstream2.str() + lastPart;

  }
  //Tornado Data
  if(selectData==5)
  {
    
    if(varid==1)
      var_name = "tornado_u_";
    else if(varid==2)
      var_name = "tornado_v_";
    else if(varid==3)
      var_name = "tornado_w_";
    

    //var_name = "tornado_mag_";  

    lastPart = ".raw";
    path = "/home/soumya/Test_DataSet/Tornado/raw_files/";

    //Construct the filename as per the value of current j
    stringstream tempstream2;

    tempstream2 <<j;

    path = path + var_name;
    fileName = path + tempstream2.str() + lastPart;
    
  }

  return fileName;                
}//Function ends

//This function computes all pair pointwise mutual information of two volumes/RVs
void computePMI_Uncert(const char *firstFile,const char *secondFile)
{
	vector <float> temp;
	vector <float> temp1;
	FILE *filept[numFiles];
	float buff=0.0;	
	float zero = 0;
	int valx=0,valy=0;
	float joint_prob_xy=0.0, prob_of_x=0.0, prob_of_y=0.0;
	N = xdim*ydim*zdim;
	float maxval,minval;
	float mm,mx;
	float count;
	vector <float> temp2;
	float weightedVal=0;
	float val=0;	

	//Resize the 4D fileStorage Vector
	rawDataStorage.resize(numFiles);  
	for(int i=0;i<rawDataStorage.size();i++)
	{
		rawDataStorage[i].resize(zdim);
		for(int j=0;j<rawDataStorage[i].size();j++)
		{
			rawDataStorage[i][j].resize(ydim);
		}
	}

	// check vti
	char *ps = (char *)firstFile+strlen(firstFile)-3;
	if (strcmp(ps, "vti")==0) {
		printf("VTI\n");
		isVTI = true;
	}

	if (!isVTI) {
		//read all files in a single large 4D vector
		filept[0] = fopen(firstFile,"rb");
		filept[1] = fopen(secondFile,"rb");

		//Now Load raw data unnormalzed
		for(int i=0;i<numFiles;i++)
		{
			for(int p=0; p<zdim; p++)
			{
				for(int q=0; q<ydim; q++)
				{
					for(int r=0; r<xdim; r++)
					{
						fread(&buff,1,sizeof(buff),filept[i]);
						rawDataStorage[i][p][q].push_back(buff);
					}
				}
			}
		}

		//close the files after reading is done
		fclose(filept[0]);
		fclose(filept[1]);
	} else
	{
	    vtkSmartPointer<vtkXMLImageDataReader> reader = vtkXMLImageDataReader::New();
	    reader->SetFileName(firstFile); //TODO
	    reader->UpdateInformation();
	    reader->Update();

	    int *ext = reader->GetOutput()->GetExtent();
	    printf("Extent: %d %d %d %d %d %d\n", ext[0], ext[1], ext[2], ext[3], ext[4], ext[5]);

	    vtkSmartPointer<vtkImageData> data = reader->GetOutput() ;

	    if (ext[1]-ext[0]+1 != xdim ||
	    	ext[3]-ext[2]+1 != ydim ||
	    	ext[5]-ext[4]+1 != zdim ) {
	    	printf("size not match\n");
	    	exit(1);
	    }

	    for (int i=0; i<numFiles; i++)
		{
	    	vtkDataArray *array = data->GetPointData()->GetArray(vtiVarname[i]);
	    	bool isVector = array->GetNumberOfComponents()>1;
			for (int z=0; z<zdim; z++)
				for (int y=0; y<ydim; y++)
					for (int x=0; x<xdim; x++)
					{
				    	if (isVector) {
				    		double value[3];
				    		array->GetTuple(x+xdim*(y+ydim*z), value);
				    		rawDataStorage[i][z][y].push_back(value[vtiTupleIdx[i]]);
				    		//printf("%lg ", value[vtiTupleIdx[i]]);
				    	} else
				    	{
				    		double value;
							array->GetTuple(x+xdim*(y+ydim*z), &value);
							rawDataStorage[i][z][y].push_back(value);
							//printf("%lg ", value);
				    	}
					}


		}

	}

	//Find Max Min for First Variable
	dataMaxFirst=dataMinFirst=rawDataStorage[0][0][0][0];
	for(int p=0; p<zdim; p++)
	{
		for(int q=0; q<ydim; q++)
		{
			for(int r=0; r<xdim; r++)
			{					
				if(dataMaxFirst<rawDataStorage[0][p][q][r])
					dataMaxFirst=rawDataStorage[0][p][q][r];

				if(dataMinFirst>rawDataStorage[0][p][q][r])
					dataMinFirst=rawDataStorage[0][p][q][r];
			}
		}		
	}

	//Find Max Min for Second Variable
	dataMaxSecond=dataMinSecond=rawDataStorage[1][0][0][0];
	for(int p=0; p<zdim; p++)
	{
		for(int q=0; q<ydim; q++)
		{
			for(int r=0; r<xdim; r++)
			{						
				if(dataMaxSecond<rawDataStorage[1][p][q][r])
					dataMaxSecond=rawDataStorage[1][p][q][r];

				if(dataMinSecond>rawDataStorage[1][p][q][r])
					dataMinSecond=rawDataStorage[1][p][q][r];
			}
		}		
	}

	for(int i=0;i<ARR_DIM;i++)
	{
		Array1[i]=0;
		Array2[i]=0;

		for(int j=0;j<ARR_DIM;j++)
		{
			ArrayComb[i][j] = 0;
		}	
	}	

	for(int k=0;k<zdim;k++)
		for(int j=0;j<ydim;j++)
			for(int i=0;i<xdim;i++)
			{
				valx = (int)(((rawDataStorage[0][k][j][i]-dataMinFirst)/(dataMaxFirst-dataMinFirst))*(ARR_DIM-1));
				valy = (int)(((rawDataStorage[1][k][j][i]-dataMinSecond)/(dataMaxSecond-dataMinSecond))*(ARR_DIM-1));

				Array1[valx]++; //nx
				Array2[valy]++; //ny
				ArrayComb[valx][valy]++; //nxy
			} 	

			//Calculatre PMI Here*************************************/		
			for(int i=0;i<ARR_DIM;i++)
			{
				for(int j=0;j<ARR_DIM;j++)
				{
					joint_prob_xy = ArrayComb[i][j] /(float) N;
					prob_of_x = Array1[i]/(float)N;
					prob_of_y = Array2[j]/(float)N;

					//if(prob_of_x*prob_of_y!=0 && joint_prob_xy!=0)
					if(joint_prob_xy>=delta)
        			{
						val = log2(joint_prob_xy/(prob_of_x*prob_of_y));

						if(val!=0)
			            temp.push_back(val); 
			            else
			            temp.push_back(0);  

			            temp1.push_back(joint_prob_xy);    									
					}
					else if(joint_prob_xy==0.0)
					{
						val=0;
						temp.push_back(0.0);
						temp1.push_back(0.0) ; 		
					}
					else
					{
						val=9999;
						temp.push_back(0.0);
						temp1.push_back(0.0) ;		
					}

				}

				PMI_Vals.push_back(temp);
				Joint_Dist_Vals.push_back(temp1);
				temp.clear();
				temp1.clear();
			}	

	//Clear the temp PMI_temp
	//rawDataStorage.clear();					
}

int main(int argc, char **argv)
{
	//argv[1] = xdim
	//argv[2] = ydim
	//argv[3] = zdim
	//argv[4] = dataset number (0=combustion, 1=Ion front, 2=Isabel 3=POP Ocean)
	//argv[5] = bin size (default 256)
	//argv[6] = timestep
	string firstFileName,secondFileName,thirdFileName,combine_name,combine_name1;
	stringstream tt;
	float val=0;
	float mean=0;
	firstFileName = argv[1];
	secondFileName = argv[2];
	xdim=atoi(argv[3]);
	ydim=atoi(argv[4]);
	zdim=atoi(argv[5]);
	selectData=atoi(argv[6]);
	ARR_DIM = atoi(argv[7]);
	timestep = atoi(argv[8]);
	vtiVarname[0] = argv[9];
	vtiVarname[1] = argv[10];

	tt<<timestep;
	Array1 = (int *)malloc(ARR_DIM*sizeof(int));
	Array2 = (int *)malloc(ARR_DIM*sizeof(int));
	temp_value = (float *)malloc(ARR_DIM*sizeof(float));
	int val1=0,val2=0;
	float pmi_value=0;
	vector <float> one;
	vector <float> two;

	I11 = (float *)malloc(ARR_DIM*sizeof(float));
	I12 = (float *)malloc(ARR_DIM*sizeof(float));
	I21 = (float *)malloc(ARR_DIM*sizeof(float));
	I22 = (float *)malloc(ARR_DIM*sizeof(float));
	I32 = (float *)malloc(ARR_DIM*sizeof(float));
	I31 = (float *)malloc(ARR_DIM*sizeof(float));

	ArrayComb = (int **)malloc(ARR_DIM*sizeof(int *));
	for(int i=0;i<ARR_DIM;i++)
		ArrayComb[i] = (int *)malloc(ARR_DIM*sizeof(int));

	CosineSimilariy = (float **)malloc(ARR_DIM*sizeof(float *));
	for(int i=0;i<ARR_DIM;i++)
		CosineSimilariy[i] = (float *)malloc(ARR_DIM*sizeof(float));


	//firstFileName = "/home/soumya/Desktop/plume/plume_zvel.raw";
	//secondFileName = "/home/soumya/Desktop/plume/gradient_vol.raw";

	//firstFileName = constructFileName(1,timestep);	
	//secondFileName = constructFileName(2,timestep);

	//firstFileName = argv[7];//constructFileName(1,timestep);	
	//secondFileName = argv[8];//constructFileName(11,timestep);

	cout<<firstFileName<<endl;
	cout<<secondFileName<<endl;
	combine_name1 = "pmi_2var_vol.raw";
	cout<<"output combined file: "<<combine_name1<<endl;   	

	//Compute the PMI: Output a 2D vector containing all pair PMI values
	computePMI_Uncert(firstFileName.c_str(),secondFileName.c_str());


	//compute_metrics();

	/*ofstream fp;
	fp.open("out.txt");
	float one_val=1;
	//float zero=0;
	for(int i=0;i<ARR_DIM;i++)
	{
		for(int j=0;j<ARR_DIM;j++)
		{
			val = PMI_Vals[i][j];
			fp<<val<<" ";
		}

		fp<<endl;
	}
	fp.close();*/

	//For printing sorting pmi values
	/*std::map<float,vec2>::iterator it;
	for(int i=0;i<ARR_DIM;i++)
	{
		for(int j=0;j<ARR_DIM;j++)
		{
			val = PMI_Vals[i][j];

			if(val!=0)
			{
				vec2 vector2 = vec2(i,j);
				pairMap.insert(pair<float,vec2>(val,vector2));
			}
		}
	}

	for (it=pairMap.begin(); it!=pairMap.end(); ++it)
	{
		cout<<it->first<<"    "<<dataMinFirst + (it->second.y/(float)ARR_DIM)*(dataMaxFirst-dataMinFirst)<<"   "
				<<dataMinSecond + (it->second.x/(float)ARR_DIM)*(dataMaxSecond-dataMinSecond)<<endl;
	}

	/*for (it=pairMap.begin(); it!=pairMap.end(); ++it)
	{
		cout<<it->first<<"    "<<(it->second.x/(float)ARR_DIM)<<"   "<<(it->second.y/(float)ARR_DIM)<<endl;
	}*/

	//Compute cosine similarity for all pairs	
	/*int i,j,k;	
	for(i=0;i<ARR_DIM;i++)
	{			
		for(j=0;j<ARR_DIM;j++)
		{
			one.clear();
			two.clear();
			for(k=0;k<ARR_DIM;k++)
			{
				one.push_back(PMI_Vals[i][k]);
				two.push_back(PMI_Vals[k][j]);
			}

			CosineSimilariy[i][j] = cosineProduct(one,two,ARR_DIM);
			//cout<<CosineSimilariy[i][j]<<endl;
		}	
		
	}*/

	//For printing sorting pmi values
	/*std::map<float,vec2>::iterator it;
	for(int i=0;i<ARR_DIM;i++)
	{
		for(int j=0;j<ARR_DIM;j++)
		{
			val = CosineSimilariy[i][j];

			if(val!=0)
			{
				vec2 vector2 = vec2(i,j);
				pairMap.insert(pair<float,vec2>(val,vector2));
			}
		}
	}

	for (it=pairMap.begin(); it!=pairMap.end(); ++it)
	{
		cout<<it->first<<"    "<<dataMinFirst + (it->second.y/(float)ARR_DIM)*(dataMaxFirst-dataMinFirst)<<"   "
				<<dataMinSecond + (it->second.x/(float)ARR_DIM)*(dataMaxSecond-dataMinSecond)<<endl;
	}

	for (it=pairMap.begin(); it!=pairMap.end(); ++it)
	{
		cout<<it->first<<"    "<<(it->second.x/(float)ARR_DIM)<<"   "<<(it->second.y/(float)ARR_DIM)<<endl;
	}*/


	//write joint pmf field
	FILE* joint_pmf_field;
	joint_pmf_field = fopen(combine_name1.c_str(),"w");
	float zero=0;
	float var_value=0;
	for(int p=0; p<zdim; p++)
	{
		for(int q=0; q<ydim; q++)
		{
			for(int r=0; r<xdim; r++)
			{           
				val1 = (rawDataStorage[0][p][q][r]-dataMinFirst)*(ARR_DIM-1)/(dataMaxFirst-dataMinFirst);
				val2 = (rawDataStorage[1][p][q][r]-dataMinSecond)*(ARR_DIM-1)/(dataMaxSecond-dataMinSecond);
				pmi_value = PMI_Vals[val1][val2];

				//pmi_value = I32[val1];

				float aa = val1/(float)(ARR_DIM-1);
				float bb = val2/(float)(ARR_DIM-1);
			
				if(pmi_value!=0)  
					fwrite(&pmi_value,sizeof(float),1,joint_pmf_field);
				else
					fwrite(&zero,sizeof(float),1,joint_pmf_field);
			}
		}   
	}

	fclose(joint_pmf_field);

	free(Array1);
	free(Array2);
	free(temp_value);
	return 0;
}
