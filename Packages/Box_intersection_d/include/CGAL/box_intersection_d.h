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
#include <cassert>

CGAL_BEGIN_NAMESPACE

template< class RandomAccessIter1, class RandomAccessIter2,
          class Callback >
void box_intersection_all_pairs_d( 
                RandomAccessIter1 p_begin, RandomAccessIter1 p_end,
                RandomAccessIter2 i_begin, RandomAccessIter2 i_end,
                Callback& callback,
                Box_intersection_d::Topology
                         topology = Box_intersection_d::CLOSED )
{
    typedef typename std::iterator_traits<RandomAccessIter1>::value_type
        Box_type;
    typedef Box_intersection_d::Box_traits_d< Box_type >
            Box_traits;
  
    box_intersection_all_pairs_custom_d( p_begin, p_end, i_begin, i_end, 
                                         callback, Box_traits(), topology );   
}

template< class RandomAccessIter1, class RandomAccessIter2,
          class Callback, class BoxTraits >
void box_intersection_all_pairs_custom_d( 
                RandomAccessIter1 p_begin, RandomAccessIter1 p_end,
                RandomAccessIter2 i_begin, RandomAccessIter2 i_end,
                Callback& callback, BoxTraits traits,
                Box_intersection_d::Topology
                         topology = Box_intersection_d::CLOSED )
{
    if (topology == Box_intersection_d::CLOSED) {
        typedef Box_intersection_d::Box_predicate_traits_d< BoxTraits, true >
                Traits;
        box_intersection_all_pairs_custom_predicates_d(
                   p_begin, p_end, i_begin, i_end, callback, Traits() );
    } else {
        typedef Box_intersection_d::Box_predicate_traits_d< BoxTraits, false >
                Traits;
        box_intersection_all_pairs_custom_predicates_d(
                   p_begin, p_end, i_begin, i_end, callback, Traits() );
    }
}


template< class RandomAccessIter1, class RandomAccessIter2,
          class Callback, class BoxPredicateTraits >
void box_intersection_all_pairs_custom_predicates_d( 
                RandomAccessIter1 p_begin, RandomAccessIter1 p_end,
                RandomAccessIter2 i_begin, RandomAccessIter2 i_end,
                Callback& callback, BoxPredicateTraits traits )
{
    Box_intersection_d::all_pairs( p_begin, p_end, i_begin, i_end, 
                                   callback, traits );   
}

template< class RandomAccessIter, class Callback >
void box_intersection_d(RandomAccessIter p_begin, RandomAccessIter p_end,
                        Callback& callback,
                        std::size_t cutoff = 10,
                        Box_intersection_d::Topology
                                    topology = Box_intersection_d::CLOSED)
{
    typedef typename std::iterator_traits<RandomAccessIter>::value_type
            Box_type;
    typedef Box_intersection_d::Box_traits_d< Box_type >
            Box_traits;

    box_intersection_custom_d(p_begin, p_end, callback,
              Box_traits(), cutoff, topology );
}


template< class RandomAccessIter1, class RandomAccessIter2, class Callback >
void box_intersection_d(RandomAccessIter1 p_begin, RandomAccessIter1 p_end,
                        RandomAccessIter2 i_begin, RandomAccessIter2 i_end,
                        Callback& callback,
                        std::size_t cutoff = 10,
                        Box_intersection_d::Topology
                                    topology = Box_intersection_d::CLOSED,
                        Box_intersection_d::Setting
                                    setting  = Box_intersection_d::BIPARTITE)
 {
    typedef typename std::iterator_traits<RandomAccessIter1>::value_type
            Box_type;
    typedef Box_intersection_d::Box_traits_d< Box_type >
            Box_traits;

    box_intersection_custom_d(p_begin, p_end, i_begin, i_end,
                       callback, Box_traits(), cutoff, topology, setting);
}

template< class RandomAccessIter, class Callback, class BoxTraits >
void box_intersection_custom_d(
     RandomAccessIter p_begin, RandomAccessIter p_end,
     Callback& callback,
     BoxTraits box_traits,
     std::size_t cutoff = 10,
     Box_intersection_d::Topology topology = Box_intersection_d::CLOSED)
{
    typedef typename std::iterator_traits<RandomAccessIter>::value_type
        Box_type;
    std::vector< Box_type > i( p_begin, p_end );
    
    box_intersection_custom_d( p_begin, p_end, i.begin(), i.end(),
      callback, box_traits, cutoff, topology, Box_intersection_d::COMPLETE );
}


template< class RandomAccessIter1, class RandomAccessIter2,
          class Callback, class BoxTraits >
void box_intersection_custom_d(
     RandomAccessIter1 p_begin, RandomAccessIter1 p_end,
     RandomAccessIter2 i_begin, RandomAccessIter2 i_end,
     Callback& callback,
     BoxTraits box_traits,
     std::size_t cutoff = 10,
     Box_intersection_d::Topology topology = Box_intersection_d::CLOSED,
     Box_intersection_d::Setting  setting  = Box_intersection_d::BIPARTITE)
 {
    if (topology == Box_intersection_d::CLOSED) {
        typedef Box_intersection_d::Box_predicate_traits_d< BoxTraits, true >
                Traits;
        box_intersection_custom_predicates_d(p_begin, p_end, i_begin, i_end,
                                        callback, Traits(), cutoff, setting);
    } else {
        typedef Box_intersection_d::Box_predicate_traits_d< BoxTraits, false >
                Traits;
        box_intersection_custom_predicates_d(p_begin, p_end, i_begin, i_end,
                                        callback, Traits(), cutoff, setting);
    }
}

template< class RandomAccessIter1, class RandomAccessIter2,
          class Callback, class BoxPredicateTraits >
void box_intersection_custom_predicates_d(
      RandomAccessIter1 p_begin, RandomAccessIter1 p_end,
      RandomAccessIter2 i_begin, RandomAccessIter2 i_end,
      Callback& callback,
      BoxPredicateTraits traits,
      std::size_t cutoff = 10,
      Box_intersection_d::Setting setting = Box_intersection_d::BIPARTITE)
{
    typedef BoxPredicateTraits Traits;
    typedef typename Traits::NT NT;
    assert( Traits::get_dim() > 0 );
    const std::size_t dim = Traits::get_dim() - 1;
    const NT inf = Box_intersection_d::box_limits<NT>::inf();
    const NT sup = Box_intersection_d::box_limits<NT>::sup();
    Box_intersection_d::segment_tree(p_begin, p_end, i_begin, i_end,
                               inf, sup, callback, traits, cutoff, dim, true);
    if(setting == Box_intersection_d::BIPARTITE)
        Box_intersection_d::segment_tree(i_begin, i_end, p_begin, p_end,
                              inf, sup, callback, traits, cutoff, dim, false);
}

CGAL_END_NAMESPACE


#endif
