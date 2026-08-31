// Copyright (c) 2026  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : François Protais

#ifndef CGAL_MESH_SMOOTHING_3_INTERNAL_LOOPS_H
#define CGAL_MESH_SMOOTHING_3_INTERNAL_LOOPS_H

#include <CGAL/license/Mesh_smoothing_3.h>

#include <CGAL/tags.h>

#include <cstddef>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>
#endif

namespace CGAL {
namespace Mesh_smoothing_3_internal {

/*
 * Reduction requirements:
 *
 *   - Reduction{} is the identity element.
 *   - reduction.join(other) combines `other` into `reduction`.
 *
 * Function must be callable as:
 *
 *   function(std::size_t index, Reduction& local_reduction);
 */
template <typename ConcurrencyTag, typename Reduction, typename Function>
Reduction reduce(std::size_t first, std::size_t last, Function&& function) {
  if constexpr(!ConcurrencyTag::is_parallel) {
    Reduction result{};

    for(std::size_t i = first; i < last; ++i)
      function(i, result);

    return result;
  } else {
#ifdef CGAL_LINKED_WITH_TBB

    return tbb::parallel_reduce(
        tbb::blocked_range<std::size_t>(first, last), Reduction{},
        [&function](const tbb::blocked_range<std::size_t>& range, Reduction local) {
          for(std::size_t i = range.begin(); i != range.end(); ++i)
            function(i, local);

          return local;
        },
        [](Reduction lhs, const Reduction& rhs) {
          lhs.join(rhs);
          return lhs;
        });

#elif defined(_OPENMP) && _OPENMP >= 201307

#pragma omp declare reduction(cgal_mesh_smoothing_join:Reduction : omp_out.join(omp_in))                               \
    initializer(omp_priv = Reduction{})

    Reduction result{};

    // using a signed type for old openMP versions
#pragma omp parallel for reduction(cgal_mesh_smoothing_join : result)
    for (std::ptrdiff_t iter_t = static_cast<std::ptrdiff_t>(first); iter_t < static_cast<std::ptrdiff_t>(last); ++iter_t) {
      function(static_cast<std::size_t>(iter_t), result);
    }
    return result;

#else

    static_assert (!std::is_convertible<ConcurrencyTag, Parallel_tag>::value,
                 "Parallel_tag is enabled but neither TBB nor OpenMP is available.");

    return Reduction{};

#endif
  }
}

template <typename ConcurrencyTag, typename Reduction, typename Function>
Reduction reduce(std::size_t size, Function&& function) {
  return reduce<ConcurrencyTag, Reduction>(0, size, std::forward<Function>(function));
}

namespace internal {
  struct Empty_reduction
  {
    void join(const Empty_reduction&) {}
  };
}

/*
 * No-reduction loop.
 *
 * This intentionally delegates to reduce() so that all backend selection
 * and OpenMP/TBB-specific code remains centralized there.
 *
 * Function must be callable as:
 *
 *   function(std::size_t index);
 */
template <typename ConcurrencyTag, typename Function>
void for_each(std::size_t first, std::size_t last, Function&& function) {
  reduce<ConcurrencyTag, internal::Empty_reduction>(first, last, [&function](std::size_t i, internal::Empty_reduction&) { function(i); });
}

template <typename ConcurrencyTag, typename Function> void for_each(std::size_t size, Function&& function) {
  for_each<ConcurrencyTag>(0, size, std::forward<Function>(function));
}

} // namespace Mesh_smoothing_3_internal
} // namespace CGAL

#endif // CGAL_MESH_SMOOTHING_3_INTERNAL_LOOPS_H
