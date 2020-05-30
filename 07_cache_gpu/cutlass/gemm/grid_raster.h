/******************************************************************************
 * Copyright (c) 2017, NVIDIA CORPORATION.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the NVIDIA CORPORATION nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 ******************************************************************************/

#pragma once

/**
 * \file
 * Abstraction for enumerating \p block_task within an input matrix
 */

#include <stdint.h>

#include "../util/util.h"


namespace cutlass {
namespace gemm {


/******************************************************************************
 * grid_raster_strategy
 ******************************************************************************/

/**
 * \brief Strategies for enumerating \p block_task within an input matrix
 */
struct grid_raster_strategy
{
    /// \brief Enumerants
    enum kind_t
    {
        /**
         * Default \p block_task assignment (currently ColumnMajor for N*,
         * RowMajor for TT, and TiledCohort for TN)
         */
        Default,

        /**
         * Column-major \p block_task assignment
         */
        ColumnMajor,

        /**
         * Row-major \p block_task assignment
         */
        RowMajor,

        /**
         * Two-level \p block_task assignment (both column-major)
         */
        TiledCohort,
    };
};



/******************************************************************************
 * grid_raster
 ******************************************************************************/

/**
 * \brief Abstraction for enumerating \p block_task within an input matrix
 *
 * NB: This generic class is not directly constructible.  Algorithm-specific
 * template specializations will provide the API functionality prescribed here.
 */
template <
    int                             BlockItemsY,    ///< Height in rows of a block-wide tile in matrix C
    int                             BlockItemsX,    ///< Width in columns of a block-wide tile in matrix C
    matrix_transform_t::kind_t      TransformA,     ///< View transform enumerant for matrix A
    matrix_transform_t::kind_t      TransformB,     ///< View transform enumerant for matrix B
    grid_raster_strategy::kind_t    RasterStrategy> ///< Strategy for enumerating \p block_task within an input matrix
struct grid_raster
{
    //-------------------------------------------------------------------------
    // Device API
    //-------------------------------------------------------------------------

    /// Thread block's base item coordinates (x, y) in matrix C
    int2 block_item_coords;

    /// Constructor
    grid_raster();

    /// Whether the thread block base coordinates are out-of-bounds for an m*n matrix C
    bool is_block_oob(int m, int n);


    //-------------------------------------------------------------------------
    // Grid launch API
    //-------------------------------------------------------------------------

    /// Compute the kernel grid extents (in thread blocks) for consuming an m*n matrix C
    static dim3 grid_dims(int m, int n);
};


/******************************************************************************
 * grid_raster (ColumnMajor specialization)
 ******************************************************************************/

/**
 * \brief Abstraction for enumerating \p block_task within an input matrix
 * (ColumnMajor specialization)
 *
 * Maps thread blocksin column-major fashion
 */
template <
    int                         BlockItemsY,          ///< Height in rows of a block-wide tile in matrix C
    int                         BlockItemsX,          ///< Width in columns of a block-wide tile in matrix C
    matrix_transform_t::kind_t  TransformA,         ///< View transform enumerant for matrix A
    matrix_transform_t::kind_t  TransformB>         ///< View transform enumerant for matrix B
struct grid_raster<
    BlockItemsY,
    BlockItemsX,
    TransformA,
    TransformB,
    grid_raster_strategy::Default>                   ///< Strategy for enumerating \p block_task within an input matrix
{
    //-------------------------------------------------------------------------
    // Device API
    //-------------------------------------------------------------------------

    /// Thread block's base item coordinates (x, y) in matrix C
    int2 block_item_coords;

    /// Constructor
    inline __device__
    grid_raster()
    {
        // blockDim.x is the fastest changing grid dim on current architectures
        block_item_coords = make_int2(
            BlockItemsX * blockIdx.y,
            BlockItemsY * blockIdx.x);
    }

    /// Whether the base \p block_item_coords are out-of-bounds for an m*n matrix C
    inline __device__
    bool is_block_oob(int m, int n)
    {
        // ColumnMajor never rasterizes fully out-of-bounds thread blocks
        return false;
    }

    //-------------------------------------------------------------------------
    // Grid launch API
    //-------------------------------------------------------------------------

    /// Compute the kernel grid extents (in thread blocks) for consuming an m*n matrix C
    inline __host__ __device__
    static dim3 grid_dims(int m, int n)
    {
        // blockDim.x is the fastest changing grid dim on current architectures
        return dim3(
            (m + BlockItemsY - 1) / BlockItemsY,
            (n + BlockItemsX - 1) / BlockItemsX);
    }
};


} // namespace gemm
} // namespace cutlass
