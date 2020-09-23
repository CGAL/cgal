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
#include <tbb/task_group.h>
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
    // Here is an illustration for n=2.
    //
    // Doing a R1-R2 intersection with a 2-split is doing 4 subpairs in parallel
    // a-c and b-d
    // r1[0][0] --a-- r1[0][1] --b-- r1[0][2]
    // r2[0][0] --c-- r2[0][1] --d-- r2[0][2]
    //
    // a-d and b-c
    // r1[1][0] --a-- r1[1][1] --b-- r1[1][2]
    // r2[1][0] --c-- r2[1][1] --d-- r2[1][2]
    //
    // Ranges must be duplicates since sorting is performed

    typedef typename std::iterator_traits<RandomAccessIter1>::value_type         val_t;
    typedef typename std::iterator_traits<RandomAccessIter1>::difference_type    diff_size;

    typedef std::vector<val_t>                                                   val_container;
    typedef typename val_container::iterator                                     It;

    static constexpr int n = 4;

    const diff_size r1s = std::distance(begin1, end1);
    const diff_size r2s = std::distance(begin2, end2);

    val_container range_1_copies, range_2_copies;
    range_1_copies.reserve(r1s * n);
    range_2_copies.reserve(r2s * n);

    const diff_size r1_step = r1s / n;
    const diff_size r2_step = r2s / n;

    for(int i=0; i<n; ++i)
    {
      range_1_copies.insert(range_1_copies.end(), begin1, end1);
      range_2_copies.insert(range_2_copies.end(), begin2, end2);
    }

    // for example for n=2, there's 'begin', 'mid', and 'end' but we leave out 'end' for convenience
    std::array<std::array<It, n>, n> range_1_iterators;
    std::array<std::array<It, n>, n> range_2_iterators;

    for(int i=0; i<n; ++i)
    {
      for(int j=0; j<n; ++j)
      {
        range_1_iterators[i][j] = range_1_copies.begin() + i * r1s + j * r1_step;
        range_2_iterators[i][j] = range_2_copies.begin() + i * r2s + j * r2_step;
      }
    }

    tbb::task_group g;

    for(int i=0; i<n; ++i)
    {
      for(int j=0; j<n; ++j)
      {
        // 'j' vs 'j+i', meaning for each value of 'i' we're matching r1[i] with a shifted version of r2[i]
        int r1_endi, r1_endj;
        int r2_endi, r2_endj;

        if(i == (n-1))
        {
          if(j == (n-1)) // both very last range and very last iterator
          {
            r1_endi = -1;
            r1_endj = -1;
          }
          else
          {
            r1_endi = i;
            r1_endj = j+1;
          }
        }
        else // i != (n-1)
        {
          if(j == (n-1)) // end of that range, but not of the full container
          {
            r1_endi = i+1;
            r1_endj = 0;
          }
          else
          {
            r1_endi = i;
            r1_endj = j+1;
          }
        }

        if(i == (n-1))
        {
          if(j+i == (n-1)) // both very last range and very last iterator
          {
            r2_endi = -1;
            r2_endj = -1;
          }
          else
          {
            r2_endi = i;
            r2_endj = (j+i+1)%n;
          }
        }
        else // i != (n-1)
        {
          if(j+i == (n-1)) // end of that range, but not of the full container
          {
            r2_endi = i+1;
            r2_endj = 0;
          }
          else
          {
            r2_endi = i;
            r2_endj = (j+i+1)%n;
          }
        }

        It r1_start = range_1_iterators[i][j];
        It r1_end = (r1_endi == -1) ? range_1_copies.end() : range_1_iterators[r1_endi][r1_endj];
        It r2_start = range_2_iterators[i][(j+i)%n];
        It r2_end = (r2_endi == -1) ? range_2_copies.end() : range_2_iterators[r2_endi][r2_endj];
        CGAL_assertion(range_1_copies.begin() <= r1_start && r1_start <= r1_end && r1_end <= range_1_copies.end());
        CGAL_assertion(range_2_copies.begin() <= r2_start && r2_start <= r2_end && r2_end <= range_2_copies.end());

        // Specify "copy by value" otherwise the values of iterators for next (i,j) iterations
        // become shared with different lambdas being run in parallel, and things go wrong
        g.run([=]{ Box_intersection_d::segment_tree( r1_start, r1_end, r2_start, r2_end,
                                                     inf, sup, callback, traits, cutoff, dim, in_order); });
      }
    }

    g.wait();
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
