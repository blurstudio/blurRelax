#pragma once
#include <cstddef>
#include <type_traits>
#include <xsimd/xsimd.hpp>
#include "fastRelax.h"

namespace relax_simd {

// Pick the SIMD batch used to hold one vertex (NUM_COMPS FLOATs).
// If the compiled-for arch has a batch exactly NUM_COMPS wide (eg. AVX/AVX2 for 4 doubles)
// use that so a whole vertex lives in one register. Otherwise fall back to the native
// batch (eg. SSE2/NEON: 2 doubles), and each vertex is split across several batches.
// This also keeps AVX-512 builds on 4-wide batches instead of 8-wide ones that would
// straddle two vertices
using sized_batch = xsimd::make_sized_batch_t<FLOAT, NUM_COMPS>;
using batch = typename std::conditional<
	std::is_void<sized_batch>::value,
	xsimd::batch<FLOAT>,
	sized_batch
>::type;

constexpr std::size_t WIDTH = batch::size; // FLOATs per batch
constexpr std::size_t PER_VERT = NUM_COMPS / WIDTH; // batches per vertex

static_assert(NUM_COMPS % WIDTH == 0, "NUM_COMPS must be a multiple of the SIMD width");

} // namespace relax_simd
