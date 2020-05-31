#pragma once
#include <stdint.h>
#include "../util/util.h"
#include "block_task.h"
#include "grid_raster.h"

namespace cutlass {
namespace gemm {

  __global__ void kernel(
                       int m,                      ///< Height in rows of op(A) and C
                       int n,                      ///< Width in columns of op(B) and C
                       int k,                      ///< Width in columns of op(A) and height in rows of op(B)
                      // epilogue_op_t op,           ///< Epilogue operation to update matrix C
                       float *d_a,               ///< Pointer to matrix A array values
                       float *d_b,               ///< Pointer to matrix B array values
                       float *d_c)               ///< Pointer to matrix C array values
{
  typedef block_task<
    16,
    16,
    4> block_task_t;

    // Declare statically-allocated shared storage
    __shared__ typename block_task_t::scratch_storage_t smem;

    // Construct and run the task
    block_task_t(
        &smem,
        d_a,
        d_b,
        d_c,
        m,
        n,
        k).run();
}


/******************************************************************************
 * Launch configuration description returned to the caller
 ******************************************************************************/

/// Return details about the launch configuration to the caller
struct launch_configuration
{
    //
    // Data members
    //

    /// cudaError_t resulting from grid launch
    cudaError_t result;

    /// Kernel grid extents in thread blocks
    dim3 grid;

    /// Thread block extents in threads
    dim3 block;

    //
    // Methods
    //

    /// Constructor
    launch_configuration():
        result(cudaSuccess),
        grid(0, 0, 0),
        block(0, 0, 0) {

    }

    /// Conversion from cudaError_t
    launch_configuration(cudaError_t result):
        result(result),
        grid(0, 0, 0),
        block(0, 0, 0) {

    }

    /// Launch configuration for Cutlass kernels
    launch_configuration(
        cudaError_t result,
        dim3 grid,
        dim3 block
    ):
        result(result),
        grid(grid),
        block(block) {

    }
};


/******************************************************************************
 * Dispatch stub
 ******************************************************************************/

/**
 * GEMM dispatch stub
 *
 * This function also serves as the autotuning entrypoint to evaluate different
 * tuning parameterizations of kernel.
 */
launch_configuration dispatch(
    int             m,                              ///< Height in rows of op(A) and C
    int             n,                              ///< Width in columns of op(B) and C
    int             k,                              ///< Width in columns of op(A) and height in rows of op(B)
    float         *d_a,                           ///< Device pointer to matrix A array values
    float         *d_b,                           ///< Device pointer to matrix B array values
    float         *d_c,                           ///< Device pointer to matrix C array values
    cudaStream_t    stream = 0,                     ///< CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool            debug_synchronous = true)       ///< Whether or not to synchronize the stream after every kernel launch
                                                    ///  to check for errors.  Also causes launch configurations to be printed
                                                    ///  to the console if DEBUG is defined.  Default is \p false.
{
  typedef grid_raster<
    64,
    64>
    grid_raster_t;
  launch_configuration config;
  config.block = dim3(64);
  int dynamic_smem_bytes = 0;
  config.grid = grid_raster_t::grid_dims(m, n);
  int sm_count;
  get_sm_count(sm_count);
  gemm::kernel
    <<< config.grid,
    config.block,
    dynamic_smem_bytes,
    stream >>>(
               m,
               n,
               k,
               d_a,
               d_b,
               d_c);
  return config;
}


} // namespace gemm
} // namespace cutlass
