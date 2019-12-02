// Copyright (c) 2004  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>
//                 Andreas Meyer <ameyer@mpi-sb.mpg.de>

#ifndef CGAL_BOX_INTERSECTION_D_H
#define CGAL_BOX_INTERSECTION_D_H

#include <CGAL/license/Box_intersection_d.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Box_intersection_d/segment_tree.h>
#include <CGAL/Box_intersection_d/Box_d.h>
#include <CGAL/Box_intersection_d/Box_with_handle_d.h>
#include <CGAL/Box_intersection_d/Box_with_info_d.h>
#include <CGAL/Box_intersection_d/Box_traits_d.h>
#include <CGAL/Box_intersection_d/box_limits.h>

#include <CGAL/use.h>
#include <CGAL/tags.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_invoke.h>
#endif

#include <iterator>
#include <vector>

////////////////////////////////////////////////////////////////////////////////////////////////
/// THE CALLBACK MUST BE THREADSAFE IF YOU ARE USING THE PARALLEL MODE
////////////////////////////////////////////////////////////////////////////////////////////////

namespace CGAL {
namespace internal {

// Generic call with custom predicate traits parameter.
template< class ConcurrencyTag,
          class RandomAccessIter1, class RandomAccessIter2,
          class Callback, class Traits >
void box_intersection_segment_tree_d(
    RandomAccessIter1 begin1, RandomAccessIter1 end1,
    RandomAccessIter2 begin2, RandomAccessIter2 end2,
    Callback callback,
    const Traits& traits,
    const std::ptrdiff_t cutoff,
    const bool in_order)
{
  typedef typename Traits::NT NT;

  CGAL_assertion(Traits::dimension() > 0);
  const int dim = Traits::dimension() - 1;

  const NT inf = Box_intersection_d::box_limits<NT>::inf();
  const NT sup = Box_intersection_d::box_limits<NT>::sup();

#ifndef CGAL_LINKED_WITH_TBB
  CGAL_static_assertion_msg (!(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                             "Parallel_tag is enabled but TBB is unavailable.");
#else // CGAL_LINKED_WITH_TBB
  if(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value)
  {
    // Below takes both box ranges and split each range in 3, and then does cross products
    // to get all combinations (9 pairs).
    //
    // 3 is chosen such that we get 9 tasks in parallel, which is close to the usual number
    // of threads (8-12) on a (current) "normal" machine. This could potentially be generalized
    // by grabbing the number of available threads (n = std::thread::hardware_concurrency()) and
    // splitting in 'i' and 'j' such that i*j=n.
    //
    // The memory footprint due to having to duplicate box ranges is empirically observed
    // to be negligible (as long as a range of pointers is passed)
    //
    typedef typename std::iterator_traits<RandomAccessIter1>::value_type   val_t;
    typedef typename std::vector< val_t>::iterator                         It;

    typename std::iterator_traits<RandomAccessIter1>::difference_type r1_third = std::distance(begin1, end1) / 3;
    typename std::iterator_traits<RandomAccessIter1>::difference_type r1_two_third = 2 * r1_third;
    typename std::iterator_traits<RandomAccessIter2>::difference_type r2_third = std::distance(begin2, end2) / 3;
    typename std::iterator_traits<RandomAccessIter2>::difference_type r2_two_third = 2 * r2_third;

    CGAL_assertion(0 <= r1_third && r1_third <= r1_two_third &&
                   (r1_two_third < std::distance(begin1, end1) || std::distance(begin1, end1) == 0));
    CGAL_assertion(0 <= r2_third && r2_third <= r2_two_third &&
                   (r2_two_third < std::distance(begin2, end2) || std::distance(begin2, end2) == 0));

    std::vector< val_t> r3( begin1, end1);
    std::vector< val_t> r4( begin2, end2);
    std::vector< val_t> r5( begin1, end1);
    std::vector< val_t> r6( begin2, end2);

    RandomAccessIter1 r1_left = begin1; std::advance(r1_left, r1_third);
    RandomAccessIter2 r2_left = begin2; std::advance(r2_left, r2_third);
    It r3_left = r3.begin(); std::advance(r3_left, r1_third);
    It r4_left = r4.begin(); std::advance(r4_left, r2_third);
    It r5_left = r5.begin(); std::advance(r5_left, r1_third);
    It r6_left = r6.begin(); std::advance(r6_left, r2_third);

    RandomAccessIter1 r1_right = begin1; std::advance(r1_right, r1_two_third);
    RandomAccessIter2 r2_right = begin2; std::advance(r2_right, r2_two_third);
    It r3_right = r3.begin(); std::advance(r3_right, r1_two_third);
    It r4_right = r4.begin(); std::advance(r4_right, r2_two_third);
    It r5_right = r5.begin(); std::advance(r5_right, r1_two_third);
    It r6_right = r6.begin(); std::advance(r6_right, r2_two_third);

    tbb::parallel_invoke([&]{ Box_intersection_d::segment_tree( begin1, r1_left, begin2, r2_left, inf, sup,
                                                                callback, traits, cutoff, dim, in_order); },
                         [&]{ Box_intersection_d::segment_tree( r1_left, r1_right, r2_left, r2_right, inf, sup,
                                                                callback, traits, cutoff, dim, in_order); },
                         [&]{ Box_intersection_d::segment_tree( r1_right, end1, r2_right, end2, inf, sup,
                                                                callback, traits, cutoff, dim, in_order); },

                         // shifted 1
                         [&]{ Box_intersection_d::segment_tree( r3.begin(), r3_left, r4_left, r4_right, inf, sup,
                                                                callback, traits, cutoff, dim, in_order); },
                         [&]{ Box_intersection_d::segment_tree( r3_left, r3_right, r4_right, r4.end(), inf, sup,
                                                                callback, traits, cutoff, dim, in_order); },
                         [&]{ Box_intersection_d::segment_tree( r3_right, r3.end(), r4.begin(), r4_left, inf, sup,
                                                                callback, traits, cutoff, dim, in_order); },

                          // shifted 2
                         [&]{ Box_intersection_d::segment_tree( r5.begin(), r5_left, r6_right, r6.end(), inf, sup,
                                                                callback, traits, cutoff, dim, in_order); },
                         [&]{ Box_intersection_d::segment_tree( r5_left, r5_right, r6.begin(), r6_left, inf, sup,
                                                                callback, traits, cutoff, dim, in_order); },
                         [&]{ Box_intersection_d::segment_tree( r5_right, r5.end(), r6_left, r6_right, inf, sup,
                                                                callback, traits, cutoff, dim, in_order); } );
  }
  else
#endif // CGAL_LINKED_WITH_TBB
  {
    Box_intersection_d::segment_tree(begin1, end1, begin2, end2, inf, sup, callback, traits, cutoff, dim, in_order);
  }
}

} // namespace internal

// Generic call with custom predicate traits parameter.
template< class ConcurrencyTag = Sequential_tag,
          class RandomAccessIter1, class RandomAccessIter2,
          class Callback, class BoxPredicateTraits >
void box_intersection_custom_predicates_d(
    RandomAccessIter1 begin1, RandomAccessIter1 end1,
    RandomAccessIter2 begin2, RandomAccessIter2 end2,
    Callback callback,
    BoxPredicateTraits traits,
    std::ptrdiff_t cutoff = 10,
    Box_intersection_d::Setting setting = Box_intersection_d::BIPARTITE)
{
  internal::box_intersection_segment_tree_d<ConcurrencyTag>(begin1, end1, begin2, end2, callback, traits, cutoff, true);
  if(setting == Box_intersection_d::BIPARTITE)
    internal::box_intersection_segment_tree_d<ConcurrencyTag>(begin2, end2, begin1, end1, callback, traits, cutoff, false);
}

// Generic call with box traits parameter.
// - make all default parameters explicit overloads (workaround)
template< class ConcurrencyTag = Sequential_tag,
          class RandomAccessIter1, class RandomAccessIter2,
          class Callback, class BoxTraits >
void box_intersection_d(
    RandomAccessIter1 begin1, RandomAccessIter1 end1,
    RandomAccessIter2 begin2, RandomAccessIter2 end2,
    Callback callback,
    BoxTraits,
    std::ptrdiff_t cutoff,
    Box_intersection_d::Topology topology,
    Box_intersection_d::Setting  setting)
{
  if (topology == Box_intersection_d::CLOSED) {
    typedef Box_intersection_d::Predicate_traits_d<BoxTraits,true> Traits;
    box_intersection_custom_predicates_d<ConcurrencyTag>(begin1, end1, begin2, end2,
                                                         callback, Traits(), cutoff, setting);
  } else {
    typedef Box_intersection_d::Predicate_traits_d<BoxTraits,false> Traits;
    box_intersection_custom_predicates_d<ConcurrencyTag>(begin1, end1, begin2, end2,
                                                         callback, Traits(), cutoff, setting);
  }
}

template< class ConcurrencyTag = Sequential_tag,
          class RandomAccessIter1, class RandomAccessIter2,
          class Callback, class BoxTraits >
void box_intersection_d(
    RandomAccessIter1 begin1, RandomAccessIter1 end1,
    RandomAccessIter2 begin2, RandomAccessIter2 end2,
    Callback callback, BoxTraits box_traits, std::ptrdiff_t cutoff,
    Box_intersection_d::Topology topology)
{
    box_intersection_d<ConcurrencyTag>( begin1, end1, begin2, end2, callback, box_traits,
                                        cutoff, topology, Box_intersection_d::BIPARTITE);
}
template< class ConcurrencyTag = Sequential_tag,
          class RandomAccessIter1, class RandomAccessIter2,
          class Callback, class BoxTraits >
void box_intersection_d(
    RandomAccessIter1 begin1, RandomAccessIter1 end1,
    RandomAccessIter2 begin2, RandomAccessIter2 end2,
    Callback callback, BoxTraits box_traits, std::ptrdiff_t cutoff)
{
  box_intersection_d<ConcurrencyTag>( begin1, end1, begin2, end2, callback, box_traits,
                                      cutoff, Box_intersection_d::CLOSED,
                                      Box_intersection_d::BIPARTITE);
}

template< class ConcurrencyTag = Sequential_tag,
          class RandomAccessIter1, class RandomAccessIter2,
          class Callback, class BoxTraits >
void box_intersection_d(
    RandomAccessIter1 begin1, RandomAccessIter1 end1,
    RandomAccessIter2 begin2, RandomAccessIter2 end2,
    Callback callback, BoxTraits box_traits)
{
  box_intersection_d<ConcurrencyTag>( begin1, end1, begin2, end2, callback, box_traits,
                                      10, Box_intersection_d::CLOSED,
                                      Box_intersection_d::BIPARTITE);
}

// Specialized call with default box traits.
// - make all default parameters explicit overloads (workaround)
template< class ConcurrencyTag = Sequential_tag,
          class RandomAccessIter1, class RandomAccessIter2, class Callback >
void box_intersection_d(
    RandomAccessIter1 begin1, RandomAccessIter1 end1,
    RandomAccessIter2 begin2, RandomAccessIter2 end2,
    Callback callback, std::ptrdiff_t cutoff,
    Box_intersection_d::Topology topology,
    Box_intersection_d::Setting  setting)
{
  typedef typename std::iterator_traits<RandomAccessIter1>::value_type val_t;
  typedef Box_intersection_d::Box_traits_d< val_t>                     Box_traits;

  box_intersection_d<ConcurrencyTag>( begin1, end1, begin2, end2, callback, Box_traits(),
                                      cutoff, topology, setting);
}

template< class ConcurrencyTag = Sequential_tag,
          class RandomAccessIter1, class RandomAccessIter2, class Callback >
void box_intersection_d(
    RandomAccessIter1 begin1, RandomAccessIter1 end1,
    RandomAccessIter2 begin2, RandomAccessIter2 end2,
    Callback callback, std::ptrdiff_t cutoff,
    Box_intersection_d::Topology topology)
{
  typedef typename std::iterator_traits<RandomAccessIter1>::value_type val_t;
  typedef Box_intersection_d::Box_traits_d< val_t>                     Box_traits;

  box_intersection_d<ConcurrencyTag>( begin1, end1, begin2, end2, callback, Box_traits(),
                                      cutoff, topology, Box_intersection_d::BIPARTITE);
}
template< class ConcurrencyTag = Sequential_tag,
          class RandomAccessIter1, class RandomAccessIter2, class Callback >
void box_intersection_d(
    RandomAccessIter1 begin1, RandomAccessIter1 end1,
    RandomAccessIter2 begin2, RandomAccessIter2 end2,
    Callback callback, std::ptrdiff_t cutoff)
{
  typedef typename std::iterator_traits<RandomAccessIter1>::value_type val_t;
  typedef Box_intersection_d::Box_traits_d< val_t>  Box_traits;
  box_intersection_d<ConcurrencyTag>( begin1, end1, begin2, end2, callback, Box_traits(),
                                      cutoff, Box_intersection_d::CLOSED,
                                      Box_intersection_d::BIPARTITE);
}
template< class ConcurrencyTag = Sequential_tag,
          class RandomAccessIter1, class RandomAccessIter2, class Callback >
void box_intersection_d(
    RandomAccessIter1 begin1, RandomAccessIter1 end1,
    RandomAccessIter2 begin2, RandomAccessIter2 end2,
    Callback callback)
{
  typedef typename std::iterator_traits<RandomAccessIter1>::value_type val_t;
  typedef Box_intersection_d::Box_traits_d< val_t>  Box_traits;
  box_intersection_d<ConcurrencyTag>( begin1, end1, begin2, end2, callback, Box_traits(),
                                      10, Box_intersection_d::CLOSED,
                                      Box_intersection_d::BIPARTITE);
}

// Generic call with box traits parameter, specialized for self-intersection.
// - make all default parameters explicit overloads (workaround)
template< class ConcurrencyTag = Sequential_tag,
          class RandomAccessIter, class Callback, class BoxTraits >
void box_self_intersection_d(
    RandomAccessIter begin, RandomAccessIter end,
    Callback callback,
    BoxTraits box_traits,
    std::ptrdiff_t cutoff,
    Box_intersection_d::Topology topology)
{
  // Copying rather than calling 'box_intersection_d(begin, end, begin, end, ...'
  // is necessary because the 'std::partition' and range splits on the first range
  // would be messed up by sorts on the second range otherwise.
  typedef typename std::iterator_traits<RandomAccessIter>::value_type val_t;
  std::vector< val_t> i( begin, end);

  box_intersection_d<ConcurrencyTag>( begin, end, i.begin(), i.end(),
                                      callback, box_traits, cutoff, topology,
                                      Box_intersection_d::COMPLETE);
}

template< class ConcurrencyTag = Sequential_tag,
          class RandomAccessIter, class Callback, class BoxTraits >
void box_self_intersection_d(
    RandomAccessIter begin, RandomAccessIter end,
    Callback callback,
    BoxTraits box_traits,
    std::ptrdiff_t cutoff)
{
    return box_self_intersection_d<ConcurrencyTag>(begin, end, callback, box_traits, cutoff,
                                                   Box_intersection_d::CLOSED);
}

template< class ConcurrencyTag = Sequential_tag,
          class RandomAccessIter, class Callback, class BoxTraits >
void box_self_intersection_d(
    RandomAccessIter begin, RandomAccessIter end,
    Callback callback,
    BoxTraits box_traits)
{
  return box_self_intersection_d<ConcurrencyTag>(begin, end, callback, box_traits, 10);
}

// Specialized call with default box traits, specialized for self-intersection.
// - make all default parameters explicit overloads (workaround)
template< class ConcurrencyTag = Sequential_tag, class RandomAccessIter, class Callback >
void box_self_intersection_d(
    RandomAccessIter begin, RandomAccessIter end,
    Callback callback)
{
    typedef typename std::iterator_traits<RandomAccessIter>::value_type val_t;
    typedef Box_intersection_d::Box_traits_d< val_t>  Box_traits;
    box_self_intersection_d<ConcurrencyTag>(begin, end, callback, Box_traits());
}

template< class ConcurrencyTag = Sequential_tag, class RandomAccessIter, class Callback >
void box_self_intersection_d(
    RandomAccessIter begin, RandomAccessIter end,
    Callback callback,
    std::ptrdiff_t cutoff)
{
    typedef typename std::iterator_traits<RandomAccessIter>::value_type val_t;
    typedef Box_intersection_d::Box_traits_d< val_t>  Box_traits;
    box_self_intersection_d<ConcurrencyTag>(begin, end, callback, Box_traits(), cutoff);
}

  template< class ConcurrencyTag = Sequential_tag, class RandomAccessIter, class Callback >
void box_self_intersection_d(
    RandomAccessIter begin, RandomAccessIter end,
    Callback callback,
    std::ptrdiff_t cutoff,
    Box_intersection_d::Topology topology)
{
    typedef typename std::iterator_traits<RandomAccessIter>::value_type val_t;
    typedef Box_intersection_d::Box_traits_d< val_t>  Box_traits;
    box_self_intersection_d<ConcurrencyTag>(begin, end, callback,
                                            Box_traits(), cutoff, topology );
}

// Generic call for trivial all-pairs algorithm with box traits parameter.
// - make all default parameters explicit overloads (workaround)
template< class ForwardIter1, class ForwardIter2,
          class Callback, class BoxTraits >
void box_intersection_all_pairs_d(
    ForwardIter1 begin1, ForwardIter1 end1,
    ForwardIter2 begin2, ForwardIter2 end2,
    Callback callback, BoxTraits)
{
    typedef Box_intersection_d::Predicate_traits_d<BoxTraits,true> Traits;
    Box_intersection_d::all_pairs( begin1, end1, begin2, end2,
                                   callback, Traits());
}

template< class ForwardIter1, class ForwardIter2,
          class Callback, class BoxTraits >
void box_intersection_all_pairs_d(
    ForwardIter1 begin1, ForwardIter1 end1,
    ForwardIter2 begin2, ForwardIter2 end2,
    Callback callback, BoxTraits,
    Box_intersection_d::Topology topology,
    Box_intersection_d::Setting setting)
{
    bool complete_case = (setting != Box_intersection_d::BIPARTITE);
    if (topology == Box_intersection_d::CLOSED) {
        typedef Box_intersection_d::Predicate_traits_d<BoxTraits,true> Traits;
        Box_intersection_d::all_pairs( begin1, end1, begin2, end2,
                                       callback, Traits(), complete_case);
    } else {
        typedef Box_intersection_d::Predicate_traits_d<BoxTraits,false> Traits;
        Box_intersection_d::all_pairs( begin1, end1, begin2, end2,
                                       callback, Traits(), complete_case);
    }
}

template< class ForwardIter1, class ForwardIter2,
          class Callback, class BoxTraits >
void box_intersection_all_pairs_d(
    ForwardIter1 begin1, ForwardIter1 end1,
    ForwardIter2 begin2, ForwardIter2 end2,
    Callback callback, BoxTraits traits,
    Box_intersection_d::Topology topology)
{
    box_intersection_all_pairs_d( begin1, end1, begin2, end2, callback, traits,
                                  topology, Box_intersection_d::BIPARTITE);
}

// Specialized call for trivial all-pairs algorithm with default box traits.
// - make all default parameters explicit overloads (workaround)
template< class ForwardIter1, class ForwardIter2, class Callback >
void box_intersection_all_pairs_d(
    ForwardIter1 begin1, ForwardIter1 end1,
    ForwardIter2 begin2, ForwardIter2 end2,
    Callback callback)
{
    typedef typename std::iterator_traits<ForwardIter1>::value_type val_t;
    typedef Box_intersection_d::Box_traits_d< val_t>  Box_traits;
    box_intersection_all_pairs_d( begin1, end1, begin2, end2,
                                  callback, Box_traits(),
                                  Box_intersection_d::CLOSED );
}

template< class ForwardIter1, class ForwardIter2, class Callback >
void box_intersection_all_pairs_d(
    ForwardIter1 begin1, ForwardIter1 end1,
    ForwardIter2 begin2, ForwardIter2 end2,
    Callback callback,
    Box_intersection_d::Topology topology)
{
    typedef typename std::iterator_traits<ForwardIter1>::value_type val_t;
    typedef Box_intersection_d::Box_traits_d< val_t>  Box_traits;
    box_intersection_all_pairs_d( begin1, end1, begin2, end2,
                                  callback, Box_traits(), topology);
}

template< class ForwardIter1, class ForwardIter2, class Callback >
void box_intersection_all_pairs_d(
    ForwardIter1 begin1, ForwardIter1 end1,
    ForwardIter2 begin2, ForwardIter2 end2,
    Callback callback,
    Box_intersection_d::Topology topology,
    Box_intersection_d::Setting  setting)
{
    typedef typename std::iterator_traits<ForwardIter1>::value_type val_t;
    typedef Box_intersection_d::Box_traits_d< val_t>  Box_traits;
    box_intersection_all_pairs_d( begin1, end1, begin2, end2,
                                  callback, Box_traits(), topology, setting);
}

// Generic call for trivial all-pairs algorithm with box traits parameter
// specialized for self-intersection test.
// - make all default parameters explicit overloads (workaround)
template< class ForwardIter, class Callback, class BoxTraits >
void box_self_intersection_all_pairs_d(
  ForwardIter begin1, ForwardIter end1, Callback callback, BoxTraits /* traits */)
{
    typedef Box_intersection_d::Predicate_traits_d<BoxTraits,true> Traits;
    Box_intersection_d::all_pairs( begin1, end1, callback, Traits());
}

template< class ForwardIter, class Callback, class BoxTraits >
void box_self_intersection_all_pairs_d(
    ForwardIter begin1, ForwardIter end1, Callback callback, BoxTraits,
    Box_intersection_d::Topology topology)
{
    if (topology == Box_intersection_d::CLOSED) {
        typedef Box_intersection_d::Predicate_traits_d<BoxTraits,true> Traits;
        Box_intersection_d::all_pairs( begin1, end1, callback, Traits());
    } else {
        typedef Box_intersection_d::Predicate_traits_d<BoxTraits,false> Traits;
        Box_intersection_d::all_pairs( begin1, end1, callback, Traits());
    }
}

// Specialized call for trivial all-pairs algorithm with default box traits.
// specialized for self-intersection test.
// - make all default parameters explicit overloads (workaround)
template< class ForwardIter, class Callback >
void box_self_intersection_all_pairs_d(
    ForwardIter begin1, ForwardIter end1, Callback callback)
{
    typedef typename std::iterator_traits<ForwardIter>::value_type val_t;
    typedef Box_intersection_d::Box_traits_d< val_t>  Box_traits;
    box_self_intersection_all_pairs_d( begin1, end1, callback, Box_traits(),
                                       Box_intersection_d::CLOSED );
}

template< class ForwardIter, class Callback >
void box_self_intersection_all_pairs_d(
    ForwardIter begin1, ForwardIter end1, Callback callback,
    Box_intersection_d::Topology topology)
{
    typedef typename std::iterator_traits<ForwardIter>::value_type val_t;
    typedef Box_intersection_d::Box_traits_d< val_t>  Box_traits;
    box_self_intersection_all_pairs_d( begin1, end1, callback, Box_traits(),
                                       topology);
}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
