#pragma once
#include "sar_bp.cu"
/*
__device__
void linear_interpolation( const float3 * a_signal,
   			         float3 * a_interp_grid,
			   const int      a_length_signal,
			   const int      a_length_interp, 
			   const float    a_conv_factor,
			   const float    a_offset,
			   const int      i_signal,
			   const int      i_interp);
*/
__device__ __inline__ 
void calc_distance( const float3 * a_pos_bp,
		    const float3   a_pos_pc,
		          float3 * a_result,
	            const int      a_size,
		    const int      i_thread);

template< unsigned int num_pixels, unsigned int num_ft, unsigned int num_pulses>
__global__
void sar_backprojection( const float3 * a_rp_data,
                         const float3 * a_pos_tx,
			 const float3 * a_pos_rx,
			 const float3 * a_pos_bp,
			       float3 * a_img_bp,
			 const float  * a_r_min,
			 const float  * a_dr,
			 const float  * a_fc,
			 const int      a_num_pulses,
			 const int      a_num_ft,
			 const int      a_num_pixels);
