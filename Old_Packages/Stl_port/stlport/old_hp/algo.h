/*
 *
 * Copyright (c) 1994
 * Hewlett-Packard Company
 *
 * Copyright (c) 1996,1997
 * Silicon Graphics Computer Systems, Inc.
 *
 * Copyright (c) 1997
 * Moscow Center for SPARC Technology
 *
 * Copyright (c) 1999 
 * Boris Fomitchev
 *
 * This material is provided "as is", with absolutely no warranty expressed
 * or implied. Any use is at your own risk.
 *
 * Permission to use or copy this software for any purpose is hereby granted 
 * without fee, provided the above notices are retained on all copies.
 * Permission to modify the code and to distribute modified code is granted,
 * provided the above notices are retained, and a notice that the code was
 * modified is included with the above copyright notice.
 *
 */

#ifndef __SGI_STL_ALGO_H
#define __SGI_STL_ALGO_H

# ifndef __SGI_STL_ALGOBASE_H
#  include <algobase.h>
# endif

# ifndef __SGI_STL_TEMPBUF_H
#  include <tempbuf.h>
# endif

# ifndef __SGI_STL_INTERNAL_HEAP_H
#  include <stl_heap.h>
# endif

# ifndef __SGI_STL_ITERATOR_H
#  include <iterator.h>
# endif

# ifndef __SGI_STL_INTERNAL_ALGO_H
#  include <stl_algo.h>
# endif

# ifndef __SGI_STL_NUMERIC_H
#  include <stl_numeric.h>
# endif

#ifdef __STL_USE_NAMESPACES

# ifdef __STL_BROKEN_USING_DIRECTIVE
using namespace __STLPORT_STD;
# else
// Names from <stl_algo.h>
using __STLPORT_STD::for_each; 
using __STLPORT_STD::find; 
using __STLPORT_STD::find_if; 
using __STLPORT_STD::adjacent_find; 
using __STLPORT_STD::count; 
using __STLPORT_STD::count_if; 
using __STLPORT_STD::search; 
using __STLPORT_STD::search_n; 
using __STLPORT_STD::swap_ranges; 
using __STLPORT_STD::transform; 
using __STLPORT_STD::replace; 
using __STLPORT_STD::replace_if; 
using __STLPORT_STD::replace_copy; 
using __STLPORT_STD::replace_copy_if; 
using __STLPORT_STD::generate; 
using __STLPORT_STD::generate_n; 
using __STLPORT_STD::remove; 
using __STLPORT_STD::remove_if; 
using __STLPORT_STD::remove_copy; 
using __STLPORT_STD::remove_copy_if; 
using __STLPORT_STD::unique; 
using __STLPORT_STD::unique_copy; 
using __STLPORT_STD::reverse; 
using __STLPORT_STD::reverse_copy; 
using __STLPORT_STD::rotate; 
using __STLPORT_STD::rotate_copy; 
using __STLPORT_STD::random_shuffle; 
using __STLPORT_STD::random_sample; 
using __STLPORT_STD::random_sample_n; 
using __STLPORT_STD::partition; 
using __STLPORT_STD::stable_partition; 
using __STLPORT_STD::sort; 
using __STLPORT_STD::stable_sort; 
using __STLPORT_STD::partial_sort; 
using __STLPORT_STD::partial_sort_copy; 
using __STLPORT_STD::nth_element; 
using __STLPORT_STD::lower_bound; 
using __STLPORT_STD::upper_bound; 
using __STLPORT_STD::equal_range; 
using __STLPORT_STD::binary_search; 
using __STLPORT_STD::merge; 
using __STLPORT_STD::inplace_merge; 
using __STLPORT_STD::includes; 
using __STLPORT_STD::set_union; 
using __STLPORT_STD::set_intersection; 
using __STLPORT_STD::set_difference; 
using __STLPORT_STD::set_symmetric_difference; 
using __STLPORT_STD::min_element; 
using __STLPORT_STD::max_element; 
using __STLPORT_STD::next_permutation; 
using __STLPORT_STD::prev_permutation; 
using __STLPORT_STD::find_first_of; 
using __STLPORT_STD::find_end; 
using __STLPORT_STD::is_sorted; 
using __STLPORT_STD::is_heap; 

// Names from stl_heap.h
using __STLPORT_STD::push_heap;
using __STLPORT_STD::pop_heap;
using __STLPORT_STD::make_heap;
using __STLPORT_STD::sort_heap;

// Names from <stl_numeric.h>
using __STLPORT_STD::accumulate; 
using __STLPORT_STD::inner_product; 
using __STLPORT_STD::partial_sum; 
using __STLPORT_STD::adjacent_difference; 
using __STLPORT_STD::power; 
using __STLPORT_STD::iota; 

# endif /* __STL_BROKEN_USING_DIRECTIVE */
#endif /* __STL_USE_NAMESPACES */

#endif /* __SGI_STL_ALGO_H */

// Local Variables:
// mode:C++
// End:
