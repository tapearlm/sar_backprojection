#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <cstdio>
#include <cstdlib>
#include <vector>

std::vector<float>
read_dr_vec()
{

  std::cout << "Reading dr_vec" << std::endl;
  char* buffer;
  long lSize;
  size_t result;
  FILE* fp = std::fopen("dr_vec.bin", "r");
  if(!fp) {
    std::perror("File opening failed");
  }

  // obtain file size:
  fseek (fp , 0 , SEEK_END);
  lSize = ftell (fp);
  rewind (fp);

  // allocate memory to contain the whole file:
  buffer = (char*) malloc (sizeof(char)*lSize);
  if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

  // copy the file into the buffer:
  result = fread (buffer,1,lSize,fp);
  if (result != lSize) {fputs ("Reading error",stderr); exit (3);}

  double * buffer_d = reinterpret_cast<double*>(buffer);
  
  size_t num_elems = lSize / sizeof(double);

  std::vector<float> dr_vec(buffer_d, buffer_d + num_elems);


  // terminate
  fclose (fp);
  free(buffer);

  return dr_vec;
}

std::vector<float>
read_fc_vec()
{

  std::cout << "Reading fc_vec" << std::endl;
  char* buffer;
  long lSize;
  size_t result;
  FILE* fp = std::fopen("fc_vec.bin", "r");
  if(!fp) {
    std::perror("File opening failed");
  }

  // obtain file size:
  fseek (fp , 0 , SEEK_END);
  lSize = ftell (fp);
  rewind (fp);

  // allocate memory to contain the whole file:
  buffer = (char*) malloc (sizeof(char)*lSize);
  if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

  // copy the file into the buffer:
  result = fread (buffer,1,lSize,fp);
  if (result != lSize) {fputs ("Reading error",stderr); exit (3);}

  double * buffer_d = reinterpret_cast<double*>(buffer);

  size_t num_elems = lSize / sizeof(double);
  std::vector<float> fc_vec(buffer_d, buffer_d + num_elems);

  // terminate
  fclose (fp);
  free(buffer);
  return fc_vec;

}

std::vector<float>
read_rmin()
{

  std::cout << "Reading rmin" << std::endl;
  char* buffer;
  long lSize;
  size_t result;
  FILE* fp = std::fopen("rmin.bin", "r");
  if(!fp) {
    std::perror("File opening failed");
  }

  // obtain file size:
  fseek (fp , 0 , SEEK_END);
  lSize = ftell (fp);
  rewind (fp);

  // allocate memory to contain the whole file:
  buffer = (char*) malloc (sizeof(char)*lSize);
  if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

  // copy the file into the buffer:
  result = fread (buffer,1,lSize,fp);
  if (result != lSize) {fputs ("Reading error",stderr); exit (3);}

  double * buffer_d = reinterpret_cast<double*>(buffer);
  size_t num_elems = lSize / sizeof(double);

  std::vector<float> rmin_vec(buffer_d, buffer_d + num_elems);
  // terminate
  fclose (fp);
  free(buffer);
  return rmin_vec;
}


std::vector<float3>
read_pos_tx()
{

  std::cout << "Reading pos_tx" << std::endl;
  char* buffer;
  long lSize;
  size_t result;
  std::vector<float3> pos_tx_vec;
  FILE* fp = std::fopen("pos_tx.bin", "r");
  if(!fp) {
    std::perror("File opening failed");
  }

  // obtain file size:
  fseek (fp , 0 , SEEK_END);
  lSize = ftell (fp);
  rewind (fp);

  // allocate memory to contain the whole file:
  buffer = (char*) malloc (sizeof(char)*lSize);
  if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

  // copy the file into the buffer:
  result = fread (buffer,1,lSize,fp);
  if (result != lSize) {fputs ("Reading error",stderr); exit (3);}
  size_t num_elems = lSize / sizeof(double);
  double * buffer_d = reinterpret_cast<double*>(buffer);

  // pos_tx are x,y,z coordinates, store in the format
  // [x0,y0,z0,x1,y1,z1,x2,y2,z2....]. Convert to float3 values instead

  pos_tx_vec.resize(num_elems/3);

  for( int ii = 0; ii < num_elems/3; ii++ )
    {

      int i_buff = ii * 3;
      pos_tx_vec[ii].x = (float)buffer_d[i_buff];
      pos_tx_vec[ii].y = (float)buffer_d[i_buff+1];
      pos_tx_vec[ii].z = (float)buffer_d[i_buff+2];
    }

  // terminate
  fclose (fp);
  free(buffer);
  return pos_tx_vec;
}

std::vector<float3>
read_pos_rx()
{

  std::cout << "Reading pos_rx" << std::endl;
  char* buffer;
  long lSize;
  size_t result;
  std::vector<float3> pos_rx_vec;
  FILE* fp = std::fopen("pos_rx.bin", "r");
  if(!fp) {
    std::perror("File opening failed");
  }

  // obtain file size:
  fseek (fp , 0 , SEEK_END);
  lSize = ftell (fp);
  rewind (fp);

  // allocate memory to contain the whole file:
  buffer = (char*) malloc (sizeof(char)*lSize);
  if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

  // copy the file into the buffer:
  result = fread (buffer,1,lSize,fp);
  if (result != lSize) {fputs ("Reading error",stderr); exit (3);}
  size_t num_elems = lSize / sizeof(double);
  double * buffer_d = reinterpret_cast<double*>(buffer);

  // pos_tx are x,y,z coordinates, store in the format
  // [x0,y0,z0,x1,y1,z1,x2,y2,z2....]. Convert to float3 values instead

  pos_rx_vec.resize(num_elems/3);

  for( int ii = 0; ii < num_elems/3; ii++ )
    {

      int i_buff = ii * 3;
      pos_rx_vec[ii].x = (float)buffer_d[i_buff];
      pos_rx_vec[ii].y = (float)buffer_d[i_buff+1];
      pos_rx_vec[ii].z = (float)buffer_d[i_buff+2];
    }

  // terminate
  fclose (fp);
  free(buffer);
  return pos_rx_vec;
}

std::vector<float3>
read_pos_bp()
{

  std::cout << "Reading pos_bp" << std::endl;
  char* buffer;
  long lSize;
  size_t result;
  std::vector<float3> pos_bp_vec;
  FILE* fp = std::fopen("pos_bp.bin", "r");
  if(!fp) {
    std::perror("File opening failed");
  }

  // obtain file size:
  fseek (fp , 0 , SEEK_END);
  lSize = ftell (fp);
  rewind (fp);

  // allocate memory to contain the whole file:
  buffer = (char*) malloc (sizeof(char)*lSize);
  if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

  // copy the file into the buffer:
  result = fread (buffer,1,lSize,fp);
  if (result != lSize) {fputs ("Reading error",stderr); exit (3);}
  size_t num_elems = lSize / sizeof(double);
  double * buffer_d = reinterpret_cast<double*>(buffer);

  // pos_tx are x,y,z coordinates, store in the format
  // [x0,y0,z0,x1,y1,z1,x2,y2,z2....]. Convert to float3 values instead

  pos_bp_vec.resize(num_elems/3);

  for( int ii = 0; ii < num_elems/3; ii++ )
    {

      int i_buff = ii * 3;
      pos_bp_vec[ii].x = (float)buffer_d[i_buff];
      pos_bp_vec[ii].y = (float)buffer_d[i_buff+1];
      pos_bp_vec[ii].z = (float)buffer_d[i_buff+2];
    }

  // terminate
  fclose (fp);
  free(buffer);
  return pos_bp_vec;
}

std::vector<float3>
read_range_profile_zp()
{

  std::cout << "Reading range_profile_zp" << std::endl;
  char* buffer;
  long lSize;
  size_t result;
  std::vector<float3> rp_data_vec;
  FILE* fp = std::fopen("range_profile_zp.bin", "r");
  if(!fp) {
    std::perror("File opening failed");
  }

  // obtain file size:
  fseek (fp , 0 , SEEK_END);
  lSize = ftell (fp);
  rewind (fp);

  // allocate memory to contain the whole file:
  buffer = (char*) malloc (sizeof(char)*lSize);
  if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

  // copy the file into the buffer:
  result = fread (buffer,1,lSize,fp);
  if (result != lSize) {fputs ("Reading error",stderr); exit (3);}
  size_t num_elems = lSize / sizeof(double);
  double * buffer_d = reinterpret_cast<double*>(buffer);

  // pos_tx are x,y,z coordinates, store in the format
  // [x0,y0,z0,x1,y1,z1,x2,y2,z2....]. Convert to float3 values instead

  rp_data_vec.resize(num_elems/2);
  const int imag_jump = num_elems/2;
  for( int ii = 0; ii < num_elems/2; ii++ )
    {

      rp_data_vec[ii].y = (float)buffer_d[ii];
      rp_data_vec[ii].z = (float)buffer_d[ii+imag_jump];
    }

  // terminate
  fclose (fp);
  free(buffer);

  return rp_data_vec;
}

int main(int argc, char** argv)
{
  std::vector<float> dr_vec = read_dr_vec();
  std::vector<float> fc_vec =  read_fc_vec();
  std::vector<float> rmin_vec = read_rmin();
  std::vector<float3> pos_tx_vec = read_pos_tx();
  std::vector<float3> pos_rx_vec = read_pos_rx();
  std::vector<float3> pos_bp_vec = read_pos_bp();
  std::vector<float3> rp_data_vec = read_range_profile_zp();

  
  std::cout << "Buffer values: " << std::endl;
  for(int ii = 0; ii < 30; ii++)
    {
      std::cout << dr_vec[ii] << std::endl;
    }
}










