#include <iostream>
#include <algorithm>
#include <stdio.h>
#include "helper_math.h"

//#define NUM_PIXELS 268800
//#define NUM_FT 402
//#define NUM_PULSES 345


/*
 * calc_distance
 * Given a buffer of xyz points and a single xyz point, calculate 
 * the total distance from each point in the buffer to the single point
 * and store in the x value of the results buffer. 
 * The results buffer is type float2 to allow buffer reuse in the main function
 * 
 * Inputs
 * a_pos_bp: The positions (x, y, z) coordinates to get the distance from
 * a_pos_pc: The single position (x, y, z) to get the distance from
 * a_result: The buffer to put the distances in ( stored in .x)
 * a_size:   The size of the a_pos_bp buffer
 * i_thread: The thread index to operate on
 */
__device__ __inline__ 
void calc_distance( const float3 * a_pos_bp,
		    const float3   a_pos_pc,
		          float3 * a_result,
	            const int      a_size,
		    const int      i_thread)
{
      
  if( i_thread >= a_size)
    return;
      
  float3 dist = a_pos_bp[i_thread] - a_pos_pc;	
  dist = dist * dist;
  dist.x = dist.x + dist.y + dist.z;
  a_result[i_thread].x = sqrt(dist.x);

}

/*
 * linear_interpolation
 * https://en.wikipedia.org/wiki/Linear_interpolation
 * This kernel inputs some signal, and interpolates / extrapolates
 * that signal to a more finely sampled grid using linear interpolation
 *
 * Inputs
 * a_signal:        The signal being interpolated
 * a_interp_grid:   The grid-points to interpolate the signal on
 * a_length_signal: The length of the signal
 * a_length_interp: The length of the interp grid
 * a_conv_factor:   The ratio of signal:interp grid points
 * a_offset:        The difference between the first signal point and interp grid point
 * i_signal:        The thread that indexes into the signal
 * i_interp:        The thread that indexes into the interpolation grid
 */
/*
__device__
void linear_interpolation( const float2 * a_signal,
   			         float3 * a_interp_grid,
			   const int      a_length_signal,
			   const int      a_length_interp, 
			   const float    a_conv_factor,
			   const float    a_offset,
			   const int      i_signal,
			   const int      i_interp)
{

  *
   * To give two points of a signal (x_a, y_a) and (x_b, y_b),
   * linear interpolation estimates the value of a signal point (x, y)
   * between these two points (i.e. x_a <= x <= x_b) by the formula
   * y = y_a + (x - x_a)(y_b - y_a)/(x_b - x_a)
   * This is not the least destructive interpolation method, but 
   * provides reasonable values and can be performed very quickly
   *
  

  // If walking off the array, return
  if( i_signal >= a_length_signal )
    return;

  if( i_interp >= a_length_interp )
    return;
  
  // Find the index of signal points nearest to the
  // interp grid point we're dealing with
  int i_data_a = floorf( a_conv_factor * i_interp
			 + a_conv_factor * a_offset );


  // The interpolation grid may be bigger than the signal
  // In this case, use the first or last values of the signal
  // for points on the interpolation grid that "overflow" the signal
  // Note that for these points, we're extrapolating not interpolating
  if( i_data_a > a_length_signal - 2 )
    i_data_a = a_length_signal - 2;

  if( i_data_a < 0 )
    i_data_a = 0;

  int i_data_b = i_data_a + 1;
  
   
  // Find the slope between the two signal points
  float slope =   ( a_signal[i_data_b].y - a_signal[i_data_a].y )
    / ( a_signal[i_data_b].x - a_signal[i_data_a].x );
  
  // Find how far the interpolation point is from the first signal point
  float x_distance = ( a_interp_grid[i_interp].x - a_signal[i_data_a].x );

  // Interpolated value is value of first signal point + (slope * distance)
  // from the that point
  a_interp_grid[i_interp].y = a_signal[i_data_a].y + ( x_distance * slope );

}
*/


/*
 * sar_backprojection
 * main GPU function to execute standard sar backprojection algorithm on
 * synthetic aperature radar data consisting of N pulses (slowtime), M
 * samples (fast-time), projecting onto an image of P pixels
 *
 * Inputs:
 * a_rp_data   : pointer to an MxN complex array containing range-pulse data (m)
 *               .x is time, .y is real, .z is imaginary
 * a_pos_tx    : pointer to a 3XN buffer containing transmitter positions at each pulse (m)
 *               .x is x_postion, .y is y_postion, .z is z-position
 * a_pos_rx    : pointer to a 3xN buffer containing receiver positions at each pulse (m)
 *               .x is x_postion, .y is y_postion, .z is z-position
 * a_pos_bp    : pointer to a 3XP buffer of pixel positions in backprojected image (m)
 *               .x is x_postion, .y is y_postion, .z is z-position
 * a_img_bp    : pointer to a 1XP complex buffer to hold the backprojected image
 *               .x is time, .y is real, .z is imaginary
 * a_rmin      : pointer to a 1XN buffer of per-pulse initial range to the first bin (m)
 * a_dr        : pointer to a 1XN buffer of per-pulse range bin sample spacings (m)
 * a_fc        : pointer to a 1XN buffer of per-pulse center frequency (Hz)
 * a_num_pulses: number of pulses (N)
 * a_num_ft    : number of fast-time samples (M)
 * a_num_pixels: number of pixels in the final image (P)
*/
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
			 const int      a_num_pixels)
{

  if (threadIdx.x >= a_num_ft || blockIdx.x >= a_num_pulses)
    return;

  //Constants
  const float PI = 3.1415927;
  // speed of light = 3e8 m/s
  const float C  = 300000000;

  // Used in the interpolation
  const int conv_factor = a_num_ft / a_num_pixels;

  /*
   * In this parallelization scheme, each block operates on a separate pulse
   * i.e. each column of input radar data has one block operating on it
   * Within a column, each block loops down the rows until the end
   * This code assumes number of threads per block >= number of rows of radar data
   */
  
  int i_column = blockIdx.x;

  // Allocate two shared memory buffers to "ping-pong" inputs and
  // outputs of operations between buffers
  // For large images, shared mem reuse is necessary
  __shared__ float3 s_buff_a[ num_pixels ];
  __shared__ float3 s_buff_b[ num_pixels ];

  // load current pulse into shared memory
  __shared__ float3 s_pulse[ num_ft ];

  // pulse data won't change loop to loop
  s_pulse[ threadIdx.x ] = a_rp_data[ threadIdx.x + blockDim.x * blockIdx.x ];
  
  // Load the tx and rx positions for this pulse
  // This is inefficient. Move to const mem?
  float3 tx_pos = a_pos_tx[ i_column ];
  float3 rx_pos = a_pos_rx[ i_column ];
  float  r_min  = a_r_min[ i_column ];
  float  dr     = a_dr[ i_column ];
  float  fc     = a_fc[ i_column ];

  if( blockIdx.x == 0 && threadIdx.x == 0 )
    {
      printf("tx_pos: %g %g %g\n", tx_pos.x, tx_pos.y, tx_pos.z);
      printf("rx_pos: %g %g %g\n", rx_pos.x, rx_pos.y, rx_pos.z);
      printf("r_min: %g\n", r_min);
      printf("dr: %g\n", dr);
      printf("fc: %g\n", fc);
      printf("rp_data: %g %g %g\n", s_pulse[0].x, s_pulse[0].y, s_pulse[0].z);

    }
  /*
  // number of times the block has to loop is the total length of the buffer
  // divided by the number of threads, rounded up
  const int i_loop_max = ( a_num_pixels + blockDim.x - 1 ) / blockDim.x;

  for(int i_loop = 0; i_loop < i_loop_max; i_loop++)
    {

      int i_thread = threadIdx.x + i_loop * blockDim.x;

      // Calc range from all pixels to transmitter, store in shared mem
      calc_distance(a_pos_bp, tx_pos, s_buff_a, num_pixels, i_thread);

      // Calc range from all pixels to reciever, store in shared mem
      calc_distance(a_pos_bp, rx_pos, s_buff_b, num_pixels, i_thread);
  
      // Sum tx and rx distances to get total distance tx -> target -> rx
      // for each pixel
      s_buff_a[ i_thread ] += s_buff_b[ i_thread ];
      
      __syncthreads();

      // create interpolation grid
      s_buff_b[ i_thread ] = ( s_buff_a[ i_thread ] - 2*r_min ) / (2*dr);
      __syncthreads();

      const int offset = s_buff_b[0].x - s_pulse[0].x;
  
      // do interpolation into s_buff_b here
      linear_interpolation(s_pulse,
			   s_buff_b,
			   num_ft,
			   num_pixels,
			   conv_factor,
			   offset,
			   threadIdx.x,
			   i_thread);
  
      // modulate to center frequency
      float3 frequency_shift = make_float3(1.0, 1.0, 2*PI * (fc / C));
      s_buff_a[ i_thread ] = expf( frequency_shift * s_buff_a[ i_thread ] );
      __syncthreads();

      s_buff_b[ i_thread ] = s_buff_b[ i_thread ] * s_buff_a[ i_thread ]; 
  
      // accumulate
      a_img_bp[ i_thread ] += s_buff_b[ i_thread ];
      
    }
  */
}



