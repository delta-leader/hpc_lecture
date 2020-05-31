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
 * Abstraction for exposing architecture-specific "dot-product-accumulate"
 * ISA operations
 */

#include <stdint.h>

#include "../util/util.h"


namespace cutlass {
namespace gemm {


/******************************************************************************
 * dp_accummulate
 ******************************************************************************/


/**
 * \brief Abstraction for exposing architecture-specific "dot-product-accumulate"
 * ISA operations
 *
 * Given two K-component vectors a and b having type value_t[K] and an addend c
 * of type accum_t, the "dot-product-accumulate" of type accum_t is computed
 * as d = x[0]*y[0] + x[1]*y[1] + ...  + x[K-1]*y[K-1] + c.
 *
 * We use the notation "dpK" to connote a K-component dot-product-accumulate.
 * For example, "dp1" is a simple multiply-add.
 *
 * For given pairing of value_t and accum_t types, the corresponding
 * dp_accummulate class will:
 *
 * - Define the member-type dp_vector_t as the appropriate K-component vector
 *   type needed to leverage architecture-specific "dot-product accumulate"
 *   ISA operations.
 * - Implement the corresponding dot-product operation between two dp_vector_t
 *   inputs a and b.
 *
 */
template <
    typename value_t,       ///< Component value type
    typename accum_t>       ///< Accumulator value type
struct dp_accummulate;



/// Default "dp1" dot-product-accumulate traits specialization for value_t->accum_t
template <
    typename value_t,       ///< Component value type
    typename accum_t>       ///< Accumulator value type
struct dp_accummulate
{
    /// Single-component "dp1" dot-product vector type
    typedef value_t dp_vector_t;


    /// Compute "dp1" float->float
    inline __device__
    static void mad(
        float &d,
        const float &a,
        const float &b,
        const float &c)
    {
        asm volatile ( "fma.rn.f32 %0, %1, %2, %3;\n"
            : "=f"(d) : "f"(a), "f"(b), "f"(c));
    }
};

} // namespace gemm
} // namespace cutlass

