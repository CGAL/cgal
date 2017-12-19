// Copyright (c) 2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>

#ifndef CGAL_PARTITIONED_POLYGON_2_H
#define CGAL_PARTITIONED_POLYGON_2_H

#include <CGAL/license/Partition_2.h>


#include <list>
#include <vector>
#include <CGAL/circulator.h>

namespace CGAL {


// this will order the diagonals around a vertex from previous to 
// next, which means a CW order if the polygon vertices are in CCW order
template<class Iterator, class Traits>
class Indirect_CW_diag_compare  
{
public:
   typedef typename Traits::Point_2        Point_2;
   typedef typename Traits::Orientation_2  Orig_orientation;
  
  Indirect_CW_diag_compare(){}
   Indirect_CW_diag_compare(Point_2 vertex, Iterator prev_ref, 
                         Iterator next_ref) : 
                _orientation(Traits().orientation_2_object()),
                _vertex(vertex),
                _prev_v_ref(prev_ref)
   { 
     _vertex_orientation = _orientation(Point_2(*_prev_v_ref), Point_2(vertex), Point_2(*next_ref));
   }

   bool
   operator()(Iterator d1, Iterator d2) 
   {
     Orientation d1_orientation = _orientation(Point_2(*_prev_v_ref), Point_2(_vertex), Point_2(*d1));
     Orientation d2_orientation = _orientation(Point_2(*_prev_v_ref), Point_2(_vertex), Point_2(*d2));
     Orientation d1_to_d2 = _orientation(Point_2(*d1), Point_2(_vertex), Point_2(*d2));

      // if both diagonals are on the same side of the line from previous 
      // vertex to this vertex then d1 comes before d2 (in CW order from
      // the edge (previous, vertex)) if one makes a left turn from d1 to d2

      if (d1_orientation == d2_orientation) return (d1_to_d2 == LEFT_TURN);

      // if d1 is on the line containing the edge (previous, vertex), then
      // the vertex must be a reflex vertex (otherwise the diagonal would
      // go outside the polygon) and d1 comes first if d2 is above the line
      // containing the three points, i.e., if it turns in the same
      // direction as the original edges around vertex.

      if (d1_orientation == COLLINEAR) 
         return (d2_orientation == _vertex_orientation);

      // opposite sides of the line containing the previous edge
      // (again, vertex must be reflex) or d2 is collinear and d1 isn't.
      // In either case, d1 comes first if it is below the line containing
      // the previous edge.

      return (d1_orientation != _vertex_orientation);
   }
private:
   Orig_orientation _orientation;
   Point_2          _vertex;
   Iterator         _prev_v_ref;
   Orientation      _vertex_orientation;
};

template <class Traits_>
class Partition_vertex;

// 
// requires 
//   Traits::Polygon_2
//   Traits::Point_2
//   Traits::Left_turn_2
//   Traits::Orientation_2
//

template <class Traits_>
class Partitioned_polygon_2 : 
                            public std::vector< Partition_vertex< Traits_ > >
{
public:
   typedef Traits_                                      Traits;
   typedef Partition_vertex<Traits>                     Vertex;
   typedef typename std::vector< Vertex >::iterator   Iterator;
   typedef Circulator_from_iterator<Iterator>           Circulator;
   typedef typename Traits::Polygon_2                   Subpolygon_2;
   typedef typename Traits::Point_2                     Point_2;
   typedef typename Traits::Left_turn_2                  Left_turn_2;
   typedef std::list<Circulator>                        Diagonal_list;
   typedef typename Diagonal_list::iterator             Diagonal_iterator;


   Partitioned_polygon_2() : _left_turn(Traits().left_turn_2_object())
   { }

   template <class InputIterator>
   Partitioned_polygon_2(InputIterator first, InputIterator beyond) :  
       _left_turn(Traits().left_turn_2_object())
   {
      for (; first != beyond; first++) {
         this->push_back(Vertex(*first));
      }
   }

   void insert_diagonal(Circulator v1_ref, Circulator v2_ref)  
   {
      (*v1_ref).insert_diagonal(v2_ref);
      (*v2_ref).insert_diagonal(v1_ref);
   }

   void prune_diagonals()
   {
      Circulator first(this->begin(), this->end(), this->begin());
      Circulator c = first;

      Diagonal_iterator d;
#ifdef CGAL_PARTITIONED_POLY_DEBUG
      std::cout << "pruning diagonals ..." << std::endl;
#endif
      do {
         d = (*c).diagonals_begin();
         while (d != (*c).diagonals_end()) {
            if (!diagonal_is_necessary(c, *d))
            {
#ifdef CGAL_PARTITIONED_POLY_DEBUG
               std::cout << "   removing from " << *c << " to " << **d
                         << std::endl;
#endif
               (**d).diagonal_erase(c);
               d = (*c).diagonal_erase(d);
            }
            else
            {
              d++;
            }
         }
#ifdef CGAL_PARTITIONED_POLY_DEBUG
         (*c).print_diagonals();
#endif
         (*c).reset_current_diagonal();
      }
      while (++c != first);
   }

   // the pruning is probably no longer necessary
   template <class OutputIterator>
   OutputIterator partition(OutputIterator result, bool prune)
   {
      // walk through each vertex and sort the diagonals
      Circulator first(this->begin(), this->end());
      Circulator c = first;
      Circulator next;
      Circulator prev = c;
      prev--;
      do
      {
         next = c;
         next++;
         (*c).sort_diagonals(prev, next);
#ifdef CGAL_PARTITIONED_POLY_DEBUG
         (*c).print_diagonals();
#endif
         prev = c;
      }
      while (++c != first);

      // now remove any diagonals that do not cut a reflex angle at one end
      if (prune) prune_diagonals();

#ifdef CGAL_PARTITIONED_POLY_DEBUG
      c = first;
      do
      {
         (*c).print_diagonals();
      }
      while (++c != first);
#endif

      make_polygon(first, result);
      return result;
   }

private:
   template<class OutputIterator>
   Circulator make_polygon(Circulator start, OutputIterator& result)
   {
       Subpolygon_2 new_polygon;
       Circulator next = start;
       do
       {
          new_polygon.push_back(*next);
#ifdef CGAL_PARTITIONED_POLY_DEBUG
          std::cout << "adding vertex " << *next << std::endl;
#endif
          Circulator diag;
          if ((*next).has_unused_diagonals())
          {
             diag = (*next).current_diagonal();
#ifdef CGAL_PARTITIONED_POLY_DEBUG
             std::cout << "diagonal endpoint: " << *diag << std::endl;
#endif
             (*next).advance_diagonal();
             if (diag == start)
             {
                *result = new_polygon;
                result++;
                return next;
             }
             else
             {
                next = make_polygon(next, result);
             }
          }
          else next++;
       } while (next != start);
       *result = new_polygon;
       result++;
       return next;
       // if there are no diagonals at this vertex
       //    push on the vertex
       // else if the first diagonal closes the polygon
       //    close the polygon
       //    return the current vertex (NOT the other end of the diagonal)
       // else
       //    remove the first diagonal
       //    recur, starting a new polygon at this vertex and return the
       //      vertex where the new polygon ended
       //    continue from the last vertex of the new polygon
   }

   

   bool cuts_reflex_angle(Circulator vertex_ref, Circulator diag_endpoint)
   {
      Circulator prev = vertex_ref; prev--;
      Circulator next = vertex_ref; next++;
   
      // find diag_endpoint in vertex_ref's list of diagonals
      Diagonal_iterator d_it;
      for (d_it = (*vertex_ref).diagonals_begin();
           d_it != (*vertex_ref).diagonals_end() && diag_endpoint != *d_it;
           d_it++)
      {
         prev = *d_it;
      }
      Diagonal_iterator next_d_it = d_it;
      next_d_it++;
      if (next_d_it == (*vertex_ref).diagonals_end())
      {
         next = vertex_ref;
         next++;
      }
      else
         next = *next_d_it;
   
//      return _right_turn(*prev, *vertex_ref, *next);
      return _left_turn(Point_2(*vertex_ref), Point_2(*prev), Point_2(*next));
   }

   bool diagonal_is_necessary(Circulator diag_ref1, Circulator diag_ref2) 
   {
       return (cuts_reflex_angle(diag_ref1, diag_ref2) ||
               cuts_reflex_angle(diag_ref2, diag_ref1));
   }

   Left_turn_2 _left_turn;
};

template <class Traits_>
class Partition_vertex : public Traits_::Point_2
{
  public:
    typedef Traits_                                              Traits;
    typedef typename Traits::Point_2                             Base_point;
    typedef typename Partitioned_polygon_2< Traits >::Circulator Circulator; 
  typedef Partition_vertex<Traits>                               Self;
//
//  It might be better if this were a set that used Indirect_CW_diag_compare
//  as the Compare object, but the constructor for Indirect_CW_diag_compare
//  requires prev and next pointers, which would have to be supplied to
//  the constructor for a Partition_vertex as well, which is difficult to
//  do (perhaps impossible in general since you don't know what next and
//  previous will be until the whole polygon is constructed)
//
    typedef std::list<Circulator>                         Diagonal_list;
    typedef typename Diagonal_list::iterator              Diagonal_iterator;

  //default constructor added for EPECK
  Partition_vertex(): Base_point()
  {
    current_diag = diag_endpoint_refs.end() ;
  }

  Partition_vertex(Base_point p)
    : Base_point(p) 
  { 
    current_diag = diag_endpoint_refs.end() ; 
  }

    Partition_vertex(const Partition_vertex& other)
      : Base_point(other) 
  { 
    // No need to deep copy.
    // We initialize in order to avoid problem with g++ safe STL
    current_diag = diag_endpoint_refs.end() ; 
  }

    void insert_diagonal(Circulator v_ref) 
    {
       diag_endpoint_refs.push_back(v_ref);
    }

    Diagonal_iterator diagonal_erase(Diagonal_iterator d_ref) 
    {
       return diag_endpoint_refs.erase(d_ref);
    }

    Diagonal_iterator diagonal_erase(Circulator diag_endpoint) 
    {
       Diagonal_iterator d_it = diagonals_begin();
       for (d_it = diagonals_begin(); d_it != diagonals_end() && 
                                     *d_it != diag_endpoint; d_it++) {}
       if (d_it != diagonals_end()) return diag_endpoint_refs.erase(d_it);
       return d_it;
    }

    Diagonal_iterator diagonals_begin() 
    {
       return diag_endpoint_refs.begin();
    }

    Diagonal_iterator diagonals_end() 
    {
       return diag_endpoint_refs.end();
    }

    bool has_unused_diagonals( )  
    {
       return current_diag != diag_endpoint_refs.end();
    }

    // sort the diagonals ccw around the point they have in common
    // and remove any duplicate diagonals
    void sort_diagonals(const Circulator& prev, const Circulator& next) 
    {
      diag_endpoint_refs.sort(Indirect_CW_diag_compare<Circulator,Traits>(*this, prev, next));

       diag_endpoint_refs.unique();
       current_diag = diag_endpoint_refs.begin();
    }

    void reset_current_diagonal( ) 
    {
       current_diag = diag_endpoint_refs.begin();
    }

    Circulator current_diagonal( ) const
    {  return *current_diag; }

    void advance_diagonal() 
    {
       if (current_diag != diag_endpoint_refs.end()) 
          current_diag++;
    }
   
    void print_diagonals( ) const
    {
       std::cout << "from " << *this << std::endl;
       typename std::list<Circulator>::const_iterator it;
       for (it = diag_endpoint_refs.begin();it != diag_endpoint_refs.end();
            it++)
       {
          std::cout << " to " << **it << std::endl;
       }
    }

private:
    Diagonal_list diag_endpoint_refs;
    Diagonal_iterator current_diag;
};

}

#endif // CGAL_PARTITIONED_POLYGON_2_H
