#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "bp_data/read_data.hpp"
#include "sar_bp.hpp"

__global__ 
void GPU_unit_test_device_fcns( const float3 * a_signal,
   			               float3 * a_interp_grid,
             			 const int      a_length_signal,
			         const int      a_length_interp, 
			         const float    a_conv_factor,
			         const float    a_offset)


{

  // Run first device function, linear_interpolation

  int i_loop_max = a_length_interp / a_length_signal;                                                                                             
                                                                                                                                                                     
  for(int i_loop = 0; i_loop < i_loop_max; i_loop++)                                                                                                                 
    { 
      /*
      int i_interp = threadIdx.x + i_loop * a_length_signal; 
            linear_interpolation(a_signal,
			   a_interp_grid,
			   a_length_signal,
			   a_length_interp,
			   a_conv_factor,
			   a_offset,
			   threadIdx.x,
			   i_interp);
      */
    }

}

/*
  * unit_test_linear_interp
  * 
  * This is a unit test of the linear interpolation function
  * It fills an signal and interpolation grid buffer with data
  * and runs the linear interpolation function. Data selected covers 
  * both buffer walk-off edge cases, positive and negative slopes, and 
  * positive and negative values
  */
__host__
void CPU_unit_test_device_fcns()
{

  // Specify a signal to interpolate.
  // For ease of testing, only use three points
  size_t signal_size = 3;
  float3 signal[signal_size];
  signal[0] = make_float3(0.0, 5.0, 0.0);
  signal[1] = make_float3(1.0, 10.0, 0.0);
  signal[2] = make_float3(2.0, 6.0, 0.0);

  // Specify a grid of positions to interpolate on
  size_t interp_grid_size = 9;
  float3 interp_grid[interp_grid_size];
  interp_grid[0] = make_float3(-1.0, 0.0, 0.0);
  interp_grid[1] = make_float3(-0.5, 0.0, 0.0);
  interp_grid[2] = make_float3(0.0, 0.0, 0.0);
  interp_grid[3] = make_float3(0.5, 0.0, 0.0);
  interp_grid[4] = make_float3(1.0, 0.0, 0.0);
  interp_grid[5] = make_float3(1.5, 0.0, 0.0);
  interp_grid[6] = make_float3(2.0, 0.0, 0.0);
  interp_grid[7] = make_float3(2.5, 0.0, 0.0);
  interp_grid[8] = make_float3(3.0, 0.0, 0.0);


  for(int ii = 0; ii < interp_grid_size; ii++)
    {
      std::cout << "On CPU interp_grid.x" << interp_grid[ii].x << ".y: " << interp_grid[ii].y << std::endl;
    } 

  // Conversion factor is how many signal points per interpolation point
  // This is used to select which signal point the linear interpolator uses for a
  // interpolation grid point

  float conv_factor = float(signal_size) / float(interp_grid_size);
  float offset = interp_grid[0].x - signal[0].x;

  // Create device buffers for the data
  float3* d_signal;
  float3* d_interp_grid;
  // Allocate and copy over the signal and interpolation grid
  
  cudaMalloc(&d_signal, signal_size*sizeof(float3));
  cudaMalloc(&d_interp_grid, interp_grid_size*sizeof(float3));

  cudaMemcpy(d_signal,
	     &signal,
	     signal_size*sizeof(float3),
	     cudaMemcpyHostToDevice); 

  cudaMemcpy(d_interp_grid,
	     &interp_grid,
	     interp_grid_size*sizeof(float3),
	     cudaMemcpyHostToDevice); 

  // Launch kernel
  GPU_unit_test_device_fcns<<<1, 128>>>( d_signal,
					 d_interp_grid,
					 signal_size,
					 interp_grid_size,
					 conv_factor,
					 offset);

  // Copy back data, and print results
  cudaMemcpy(interp_grid,
	     d_interp_grid,
	     interp_grid_size*sizeof(float3),
	     cudaMemcpyDeviceToHost);

  std::cout << "Results:\n";
  for(int ii = 0; ii < interp_grid_size; ii++)
    {
      std::cout << "x: " << interp_grid[ii].x
		<< "  y: " << interp_grid[ii].y << "\n";
    }

  // cleanup
  cudaFree( d_signal );
  cudaFree( d_interp_grid);

}


__host__
bool
unit_test_sar_backprojection()
{
  //template< unsigned int num_pixels, unsigned int num_ft, unsigned int num_pulses>


  std::vector<float> dr_vec = read_dr_vec();
  float* h_dr_buff = dr_vec.data();
  float* d_dr_buff;
  std::vector<float> fc_vec =  read_fc_vec();
  float* h_fc_buff = fc_vec.data();
  float* d_fc_buff;
  std::vector<float> rmin_vec = read_rmin();
  float* h_rmin_buff = rmin_vec.data();
  float* d_rmin_buff;
  std::vector<float3> pos_tx_vec = read_pos_tx();
  float3* h_pos_tx_buff = pos_tx_vec.data();
  float3* d_pos_tx_buff;
  std::vector<float3> pos_rx_vec = read_pos_rx();
  float3* h_pos_rx_buff = pos_rx_vec.data();
  float3* d_pos_rx_buff;
  std::vector<float3> pos_bp_vec = read_pos_bp();
  float3* h_pos_bp_buff = pos_bp_vec.data();
  float3* d_pos_bp_buff;
  std::vector<float3> rp_data_vec = read_range_profile_zp();
  float3* h_rp_data_buff = rp_data_vec.data();
  float3* d_rp_data_buff;
  
  std::vector<float3> img_bp;
  img_bp.resize(pos_bp_vec.size());
  float3* h_img_bp_buff = img_bp.data();
  float3 * d_img_bp_buff;
  
  const int pos_bp_size = 268800;
  const int num_ft = 402;
  const int num_pulses = 345;
  
  
  cudaMalloc(&d_dr_buff, dr_vec.size()*sizeof(float));
  cudaMalloc(&d_fc_buff, fc_vec.size()*sizeof(float));
  cudaMalloc(&d_rmin_buff, rmin_vec.size()*sizeof(float));
  cudaMalloc(&d_pos_tx_buff, pos_tx_vec.size()*sizeof(float3));
  cudaMalloc(&d_pos_rx_buff, pos_rx_vec.size()*sizeof(float3));
  cudaMalloc(&d_pos_bp_buff, pos_bp_vec.size()*sizeof(float3));
  cudaMalloc(&d_rp_data_buff, rp_data_vec.size()*sizeof(float3));
  cudaMalloc(&d_img_bp_buff, img_bp.size()*sizeof(float3));


  cudaMemcpy(d_dr_buff,
	     &(h_dr_buff),
	     dr_vec.size()*sizeof(float),
	     cudaMemcpyHostToDevice); 

  cudaMemcpy(d_fc_buff,
	     &(h_fc_buff),
	     fc_vec.size()*sizeof(float),
	     cudaMemcpyHostToDevice); 
  cudaMemcpy(d_rmin_buff,
	     &(h_rmin_buff),
	     rmin_vec.size()*sizeof(float),
	     cudaMemcpyHostToDevice); 
  cudaMemcpy(d_pos_tx_buff,
	     &(h_pos_rx_buff),
	     pos_tx_vec.size()*sizeof(float3),
	     cudaMemcpyHostToDevice); 
  cudaMemcpy(d_pos_rx_buff,
	     &(h_pos_rx_buff),
	     pos_rx_vec.size()*sizeof(float3),
	     cudaMemcpyHostToDevice); 
  cudaMemcpy(d_pos_bp_buff,
	     &(h_pos_bp_buff),
	     pos_bp_vec.size()*sizeof(float3),
	     cudaMemcpyHostToDevice); 
  cudaMemcpy(d_rp_data_buff,
	     &(h_rp_data_buff),
	     rp_data_vec.size()*sizeof(float3),
	     cudaMemcpyHostToDevice); 
  cudaMemcpy(d_img_bp_buff,
	     &(h_img_bp_buff),
	     img_bp.size()*sizeof(float3),
	     cudaMemcpyHostToDevice); 

  
  
  sar_backprojection< pos_bp_size, num_ft, num_pulses >
                     <<<1, 512>>>
                    ( d_rp_data_buff,
		      d_pos_tx_buff,
		      d_pos_rx_buff,
		      d_pos_bp_buff,
		      d_img_bp_buff,
		      d_rmin_buff,
		      d_dr_buff,
		      d_fc_buff,
		      fc_vec.size(),
		      rp_data_vec.size()/fc_vec.size(),
		      pos_bp_vec.size());


  cudaFree(d_rp_data_buff);
  cudaFree(d_pos_tx_buff);
  cudaFree(d_pos_rx_buff);
  cudaFree(d_pos_bp_buff);
  cudaFree(d_img_bp_buff);
  cudaFree(d_rmin_buff);
  cudaFree(d_dr_buff);
  cudaFree(d_fc_buff);

  
 return true;
}


int main(int argc, char** argv)
{

  //  CPU_unit_test_device_fcns();
  unit_test_sar_backprojection();

}

