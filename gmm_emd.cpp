#include <string>
#include <cstdio>
#include <vector>
#include <cassert>
#include <cstring>
#include <cstdlib>

#include <vtkSmartPointer.h>
#include <vtkXMLMultiBlockDataReader.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkStructuredGrid.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkXMLMultiBlockDataWriter.h>

#include <boost/math/distributions/normal.hpp>
#include <jclib/statistics.h>

using namespace std;

#define vsp_new(type, x) vtkSmartPointer<type> x = vtkSmartPointer<type>::New()
#define PASSAGES 36
#define BLOCKS_PER_PSG 3  // blocks per passage
#define NBINS 256
const int xdim = 10; // per block
const int ydim = 14;
const int zdim = 11;
const int psg_size = xdim*ydim*zdim*BLOCKS_PER_PSG;

float meanAry[PASSAGES][4][psg_size];
float stdAry[PASSAGES][4][psg_size];
float weightAry[PASSAGES][4][psg_size];

void merge_arrays(vtkMultiBlockDataSet *mb)
{
  for (int i=0; i<PASSAGES; i++)
  {
    for (int b=0; b<BLOCKS_PER_PSG; b++)
    {
      vtkStructuredGrid *sgrid = vtkStructuredGrid::SafeDownCast( mb->GetBlock(b+i*BLOCKS_PER_PSG) );
      int *dim = sgrid->GetDimensions();
      assert(dim[0] == xdim && dim[1] == ydim && dim[2] == zdim);
      int size = xdim*ydim*zdim;
      memcpy( &meanAry[i][0][size*b], sgrid->GetPointData()->GetArray("mean0")->GetVoidPointer(0), size*sizeof(float) );
      memcpy( &meanAry[i][1][size*b], sgrid->GetPointData()->GetArray("mean1")->GetVoidPointer(0), size*sizeof(float) );
      memcpy( &meanAry[i][2][size*b], sgrid->GetPointData()->GetArray("mean2")->GetVoidPointer(0), size*sizeof(float) );
      memcpy( &meanAry[i][3][size*b], sgrid->GetPointData()->GetArray("mean3")->GetVoidPointer(0), size*sizeof(float) );
      memcpy( &stdAry[i][0][size*b], sgrid->GetPointData()->GetArray("stdev0")->GetVoidPointer(0), size*sizeof(float) );
      memcpy( &stdAry[i][1][size*b], sgrid->GetPointData()->GetArray("stdev1")->GetVoidPointer(0), size*sizeof(float) );
      memcpy( &stdAry[i][2][size*b], sgrid->GetPointData()->GetArray("stdev2")->GetVoidPointer(0), size*sizeof(float) );
      memcpy( &stdAry[i][3][size*b], sgrid->GetPointData()->GetArray("stdev3")->GetVoidPointer(0), size*sizeof(float) );
      memcpy( &weightAry[i][0][size*b], sgrid->GetPointData()->GetArray("weight0")->GetVoidPointer(0), size*sizeof(float) );
      memcpy( &weightAry[i][1][size*b], sgrid->GetPointData()->GetArray("weight1")->GetVoidPointer(0), size*sizeof(float) );
      memcpy( &weightAry[i][2][size*b], sgrid->GetPointData()->GetArray("weight2")->GetVoidPointer(0), size*sizeof(float) );
      memcpy( &weightAry[i][3][size*b], sgrid->GetPointData()->GetArray("weight3")->GetVoidPointer(0), size*sizeof(float) );
    }
  }
}

// gmm -> cdf
void getCdf(float mAry[4], float sAry[4], float wAry[4], float edges[NBINS], float hist[NBINS])
{
  // init
  for (int b=0; b<NBINS; b++)
    hist[b] = 0;

  for (int i=0; i<4; i++)
  {
    boost::math::normal_distribution<float> gauss(mAry[i], sAry[i]);
    for (int b=0; b<NBINS; b++)
    {
      hist[b] += boost::math::cdf(gauss, edges[b]) * wAry[i];
    }
  }
#if 0
  for (int b=0; b<NBINS; b++)
  {
    printf("%f ", hist[b]);
  }
  printf("\n");
#endif
}

int main(int argc, char **argv)
{
  printf ("gmm_emd <filename> <output_file_id>\n");
  if (argc<=2)
    return 1;
  string filename = argv[1];
  int file_id = atoi(argv[2]);

  vsp_new(vtkXMLMultiBlockDataReader, reader);
  reader->SetFileName(filename.c_str());
  reader->Update();
  vtkMultiBlockDataSet *mb = vtkMultiBlockDataSet::SafeDownCast( reader->GetOutput() );
  assert(mb);

  merge_arrays(mb);

  vector<vector<float> > anomaly_ary (PASSAGES, vector<float>(psg_size));

  float ptMeanAry[PASSAGES][4];
  float ptStdAry[PASSAGES][4];
  float ptWeightAry[PASSAGES][4];

  for (int idx=0; idx < psg_size; idx++)
  {
    float minVal = 1e+30;
    float maxVal = 0;
    for (int p=0; p<PASSAGES; p++)
    {
      for (int m=0; m<4; m++)
      {
        ptMeanAry[p][m] = meanAry[p][m][idx];
        ptStdAry[p][m] = stdAry[p][m][idx];
        ptWeightAry[p][m] = weightAry[p][m][idx];
        minVal = min(minVal, ptMeanAry[p][m]-ptStdAry[p][m]*3);
        maxVal = max(maxVal, ptMeanAry[p][m]+ptStdAry[p][m]*3);
      } // model
    } // b
    minVal = 0;
    maxVal = 2.5; //!!!!!!!!!!

    // compute emd
    float bw = (maxVal - minVal)/(NBINS-1); // bin width
    float edges[NBINS] ;
    for (int bin=0; bin<NBINS; bin++)
      edges[bin] = minVal + bw * bin;

    // create cdf matrix
    float cdfmat[PASSAGES][NBINS];
    {
      for (int p=0; p<PASSAGES; p++)
      {
        getCdf(ptMeanAry[p], ptStdAry[p], ptWeightAry[p], edges, cdfmat[p]);
      }
    }

    // get mean distribution
    float mean_cdf[NBINS];
    {
      for (int i=0; i<NBINS; i++)
      {
        mean_cdf[i] = 0;
        for (int j=0; j<PASSAGES; j++)
          mean_cdf[i] += cdfmat[j][i];
        mean_cdf[i] /= PASSAGES;
      }
    }

    // compute EMD
    {
      for (int p=0; p<PASSAGES; p++)
      {
        anomaly_ary[p][idx] = 0;
        for (int i=0; i<NBINS; i++)
        {
          anomaly_ary[p][idx] += fabs(mean_cdf[i] - cdfmat[p][i]) * bw;
        }
      }
    }

  } // idx

  for (int p=0; p<PASSAGES; p++)
  {
    char out[1024];
    sprintf(out, "emd_b%d_%d.raw", p+1, file_id);
    printf("saving: %s\n", out);
    FILE *fp = fopen(out, "wb");
    fwrite(&anomaly_ary[p][0], sizeof(float), psg_size, fp);
    fclose(fp);
  }

#if 1  // save vtk format
  vsp_new(vtkXMLMultiBlockDataWriter, writer);

  for (int p=0; p<PASSAGES; p++)
  {
    for (int b=0; b<BLOCKS_PER_PSG; b++)
    {
      int size = xdim*ydim*zdim;
      vsp_new(vtkFloatArray, outAry);
      outAry->SetNumberOfComponents(1);
      outAry->SetNumberOfTuples(size);
      outAry->SetName("Anomaly");
      memcpy(outAry->GetVoidPointer(0), &anomaly_ary[p][size*b], sizeof(float)*size);
      vtkStructuredGrid::SafeDownCast( mb->GetBlock(b+p*BLOCKS_PER_PSG) )->GetPointData()->AddArray(outAry);
    }
  }

  char out[1024];
  sprintf(out, "emd_%d.vtm", file_id);
  printf("saving: %s\n", out);
  writer->SetFileName(out);
  writer->SetInputData(mb);
  writer->SetEncodeAppendedData(0);
  writer->Write();
#endif
}
