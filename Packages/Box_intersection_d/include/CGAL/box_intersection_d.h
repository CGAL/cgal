// Copyright (c) 2004  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>
//                 Andreas Meyer <ameyer@mpi-sb.mpg.de>

#ifndef CGAL_BOX_INTERSECTION_D_H
#define CGAL_BOX_INTERSECTION_D_H

#include <CGAL/Box_intersection_d/segment_tree.h>
#include <CGAL/Box_intersection_d/box_traits.h>
#include <CGAL/Box_intersection_d/box_limits.h>

#include <vector>

CGAL_BEGIN_NAMESPACE

// Generic call with custom predicate traits parameter.
template< class RandomAccessIter1, class RandomAccessIter2,
          class Callback, class BoxPredicateTraits >
void box_intersection_custom_predicates_d(
    RandomAccessIter1 begin1, RandomAccessIter1 end1,
    RandomAccessIter2 begin2, RandomAccessIter2 end2,
    Callback& callback,
    BoxPredicateTraits traits,
    std::size_t cutoff = 10,
    Box_intersection_d::Setting setting = Box_intersection_d::BIPARTITE)
{
    typedef BoxPredicateTraits Traits;
    typedef typename Traits::NT NT;
    CGAL_assertion( Traits::get_dim() > 0 );
    const std::size_t dim = Traits::get_dim() - 1;
    const NT inf = Box_intersection_d::box_limits<NT>::inf();
    const NT sup = Box_intersection_d::box_limits<NT>::sup();
    Box_intersection_d::segment_tree(begin1, end1, begin2, end2,
                              inf, sup, callback, traits, cutoff, dim, true);
    if(setting == Box_intersection_d::BIPARTITE)
        Box_intersection_d::segment_tree(begin2, end2, begin1, end1,
                              inf, sup, callback, traits, cutoff, dim, false);
}


// Generic call with box traits parameter.
template< class RandomAccessIter1, class RandomAccessIter2,
          class Callback, class BoxTraits >
void box_intersection_d(
    RandomAccessIter1 begin1, RandomAccessIter1 end1,
    RandomAccessIter2 begin2, RandomAccessIter2 end2,
    Callback& callback,
    BoxTraits box_traits,
    std::size_t cutoff = 10,
    Box_intersection_d::Topology topology = Box_intersection_d::CLOSED,
    Box_intersection_d::Setting  setting  = Box_intersection_d::BIPARTITE)
{
    if (topology == Box_intersection_d::CLOSED) {
        typedef Box_intersection_d::Box_predicate_traits_d<BoxTraits,true> Tr;
        box_intersection_custom_predicates_d(begin1, end1, begin2, end2,
                                             callback, Tr(), cutoff, setting);
    } else {
        typedef Box_intersection_d::Box_predicate_traits_d<BoxTraits,false> Tr;
        box_intersection_custom_predicates_d(begin1, end1, begin2, end2,
                                             callback, Tr(), cutoff, setting);
    }
}

// Specialized call with default box traits.
template< class RandomAccessIter1, class RandomAccessIter2, class Callback >
void box_intersection_d(
    RandomAccessIter1 begin1, RandomAccessIter1 end1,
    RandomAccessIter2 begin2, RandomAccessIter2 end2,
    Callback& callback,
    std::size_t cutoff = 10,
    Box_intersection_d::Topology topology = Box_intersection_d::CLOSED,
    Box_intersection_d::Setting  setting  = Box_intersection_d::BIPARTITE)
{
    typedef typename std::iterator_traits<RandomAccessIter1>::value_type Box_t;
    typedef Box_intersection_d::Box_traits_d< Box_t>  Box_traits;
    box_intersection_d( begin1, end1, begin2, end2,
                        callback, Box_traits(), cutoff, topology, setting);
}


// Generic call with box traits parameter, specialized for self-intersection.
template< class RandomAccessIter, class Callback, class BoxTraits >
void box_self_intersection_d(
    RandomAccessIter begin, RandomAccessIter end,
    Callback& callback,
    BoxTraits box_traits,
    std::size_t cutoff = 10,
    Box_intersection_d::Topology topology = Box_intersection_d::CLOSED)
{
    typedef typename std::iterator_traits<RandomAccessIter>::value_type Box_t;
    std::vector< Box_t> i( begin, end);
    box_intersection_d( begin, end, i.begin(), i.end(),
        callback, box_traits, cutoff, topology, Box_intersection_d::COMPLETE);
}

// Specialized call with default box traits, specialized for self-intersection.
template< class RandomAccessIter, class Callback >
void box_self_intersection_d(
    RandomAccessIter begin, RandomAccessIter end,
    Callback& callback,
    std::size_t cutoff = 10,
    Box_intersection_d::Topology
    topology = Box_intersection_d::CLOSED)
{
    typedef typename std::iterator_traits<RandomAccessIter>::value_type Box_t;
    typedef Box_intersection_d::Box_traits_d< Box_t>  Box_traits;
    box_self_intersection_d(p_begin, p_end, callback,
                            Box_traits(), cutoff, topology );
}


// Generic call for trivial all-pairs algorithm with custom predicate traits.
template< class RandomAccessIter1, class RandomAccessIter2,
          class Callback, class BoxPredicateTraits >
void box_intersection_all_pairs_custom_predicates_d(
    RandomAccessIter1 begin1, RandomAccessIter1 end1,
    RandomAccessIter2 begin2, RandomAccessIter2 end2,
    Callback& callback, BoxPredicateTraits traits )
{
    Box_intersection_d::all_pairs(begin1, end1, begin2, end2, callback,traits);
}


// Generic call for trivial all-pairs algorithm with box traits parameter.
template< class RandomAccessIter1, class RandomAccessIter2,
          class Callback, class BoxTraits >
void box_intersection_all_pairs_custom_d( 
    RandomAccessIter1 begin1, RandomAccessIter1 end1,
    RandomAccessIter2 begin2, RandomAccessIter2 end2,
    Callback& callback, BoxTraits traits,
    Box_intersection_d::Topology topology = Box_intersection_d::CLOSED )
{
    if (topology == Box_intersection_d::CLOSED) {
        typedef Box_intersection_d::Box_predicate_traits_d<BoxTraits,true> Tr;
        box_intersection_all_pairs_custom_predicates_d(
                   begin1, end1, begin2, end2, callback, Tr());
    } else {
        typedef Box_intersection_d::Box_predicate_traits_d<BoxTraits,false> Tr;
        box_intersection_all_pairs_custom_predicates_d(
                   begin1, end1, begin2, end2, callback, Tr());
    }
}

// Specialized call for trivial all-pairs algorithm with default box traits.
template< class RandomAccessIter1, class RandomAccessIter2, class Callback >
void box_intersection_all_pairs_d( 
    RandomAccessIter1 begin1, RandomAccessIter1 end1,
    RandomAccessIter2 begin2, RandomAccessIter2 end2,
    Callback& callback,
    Box_intersection_d::Topology topology = Box_intersection_d::CLOSED )
{
    typedef typename std::iterator_traits<RandomAccessIter1>::value_type Box_t;
    typedef Box_intersection_d::Box_traits_d< Box_t>  Box_traits;
    box_intersection_all_pairs_custom_d( begin1, end1, begin2, end2, 
                                         callback, Box_traits(), topology );   
}


CGAL_END_NAMESPACE


#endif
