// ============================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision $
// release_date  : $CGAL_Date $
//
// file          : include/CGAL/Partitioned_polygon_2.h
// package       : $CGAL_Package: Partition_2 $
// maintainer    : Susan Hert <hert@mpi-sb.mpg.de>
// chapter       : Planar Polygon Partitioning
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Susan Hert <hert@mpi-sb.mpg.de>
//
// coordinator   : MPI (Susan Hert <hert@mpi-sb.mpg.de>)
//
// implementation: Polygon with diagonals used to represent a partitioning
// ============================================================================

#ifndef CGAL_PARTITIONED_POLYGON_2_H
#define CGAL_PARTITIONED_POLYGON_2_H

#include <list>
#include <vector>
#include <CGAL/circulator.h>
#include <CGAL/Turn_reverser.h>

namespace CGAL {

template <class Traits> 
class Partitioned_polygon_2;

// this will order the diagonals around a vertex from previous to 
// next, which means a CW order if the polygon vertices are in CCW order
template<class Iterator, class Traits>
class Indirect_CW_diag_compare  
{
public:
   typedef typename Traits::R              R;
   typedef typename Traits::Orientation_2  Orig_orientation;

   Indirect_CW_diag_compare(Point_2<R> vertex, Iterator prev_ref, 
                         Iterator next_ref) : 
                _orientation(Traits().orientation_2_object()),
                _vertex(vertex),
                _prev_v_ref(prev_ref)
   { 
      _vertex_orientation = _orientation(*_prev_v_ref, vertex, *next_ref);
   }

   bool
   operator()(Iterator d1, Iterator d2) 
   {
      Orientation d1_orientation = _orientation(*_prev_v_ref, _vertex, *d1);
      Orientation d2_orientation = _orientation(*_prev_v_ref, _vertex, *d2);
      Orientation d1_to_d2 = _orientation(*d1, _vertex, *d2);

      // if both diagonals are on the same side of the line from previous 
      // vertex to this vertex then d1 comes before d2 (in CW order from
      // the edge (previous, vertex)) if one makes a left turn from d1 to d2

      if (d1_orientation == d2_orientation) return (d1_to_d2 == LEFTTURN);

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
   Point_2<R>       _vertex;
   Iterator         _prev_v_ref;
   Orientation      _vertex_orientation;
};

template <class Traits>
class Partition_vertex : public Traits::Point_2
{
  public:
    typedef typename Traits::Point_2                             Point_2;
    typedef typename Partitioned_polygon_2< Traits >::iterator   I; 
    typedef Circulator_from_iterator<I>                          Circulator;
//
//  It might be better if this were a set that used Indirect_CW_diag_compare
//  as the Comapre object, but the constructor for Indirect_CW_diag_compare
//  requires prev and next pointers, which would have to be supplied to
//  the constructor for a Partition_vertex as well, which is difficult to
//  do (perhaps impossible in general since you don't know what next and
//  previous will be until the whole polygon is constructed)
//
    typedef std::list<Circulator>                         Diagonal_list;
    typedef typename Diagonal_list::iterator              Diagonal_iterator;


    Partition_vertex(Point_2 p): Point_2(p) {}

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
                                     *d_it != diag_endpoint; d_it++);
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
       diag_endpoint_refs.sort(
          Indirect_CW_diag_compare<Circulator,Traits>(*this, prev, next));
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
       for (it = diag_endpoint_refs.begin();it != diag_endpoint_refs.end();it++)
       {
          std::cout << " to " << **it << std::endl;
       }
    }

private:
    Diagonal_list diag_endpoint_refs;
    Diagonal_iterator current_diag;
};

// 
// requires 
//   Traits::Polygon_2
//   Traits::Point_2
//   Traits::Leftturn_2
//   Traits::Orientation_2
//   Traits::leftturn_2_object
//   Traits::orientation_2_object
template <class Traits>
class Partitioned_polygon_2 : public std::vector< Partition_vertex< Traits > >
{
   typedef Partitioned_polygon_2<Traits>             Self;
   typedef Partition_vertex<Traits>                  Vertex;
   typedef typename Self::iterator                   Iterator;
   typedef Circulator_from_iterator<Iterator>        Circulator;
   typedef typename Traits::Polygon_2                Subpolygon_2;
   typedef typename Traits::Point_2                  Point_2;
   typedef typename Traits::Leftturn_2               Leftturn_2;
   typedef Turn_reverser<Point_2, Leftturn_2>        Rightturn_2;

public:
   Partitioned_polygon_2() : 
       _rightturn(Rightturn_2(Traits().leftturn_2_object()))
   { }

   template <class InputIterator>
   Partitioned_polygon_2(InputIterator first, InputIterator beyond) :  
       _rightturn(Rightturn_2(Traits().leftturn_2_object()))
   {
      for (; first != beyond; first++) {
         push_back(Vertex(*first));
      }
   }

   void insert_diagonal(Circulator v1_ref, Circulator v2_ref)  
   {
      (*v1_ref).insert_diagonal(v2_ref);
      (*v2_ref).insert_diagonal(v1_ref);
   }

   void prune_diagonals();

   // the pruning is probably no longer necessary
   template <class OutputIterator>
   OutputIterator partition(OutputIterator result, bool prune);

private:
   template<class OutputIterator>
   Circulator make_polygon(Circulator start, OutputIterator& result);

   bool cuts_reflex_angle(Circulator vertex_ref, Circulator diag_endpoint);

   bool diagonal_is_necessary(Circulator diag_ref1, Circulator diag_ref2) 
   {
       return (cuts_reflex_angle(diag_ref1, diag_ref2) ||
               cuts_reflex_angle(diag_ref2, diag_ref1));
   }

   Rightturn_2 _rightturn;
};

}


#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/Partitioned_polygon_2.C>
#endif // CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION

#endif // CGAL_PARTITIONED_POLYGON_2_H
