// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Constrained_Delaunay_triangulation_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$

// author(s)     : Mariette Yvinec, Jean Daniel Boissonnat
//
// coordinator   : Mariette Yvinec  < Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_WI_2_H
#define CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_WI_2_H

#include <math.h>
#include <pair.h>
#include <list.h>
#include <vector.h>
#include <map.h> 
#include <set.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Constrained_triangulation_wi_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

CGAL_BEGIN_NAMESPACE
template < class  Ctrwi >
class Constrained_Delaunay_triangulation_wi_2
  : public Constrained_Delaunay_triangulation_2<Ctrwi>
{
public:
  typedef Ctrwi                                           C_triangulation_wi;
  typedef Constrained_Delaunay_triangulation_2<Ctrwi>     CD_triangulation;
  typedef Constrained_Delaunay_triangulation_wi_2<Ctrwi>  CD_triangulation_wi;
  typedef typename C_triangulation_wi::Geom_traits        Geom_traits;
  typedef typename C_triangulation_wi::Constraint         Constraint;
  typedef typename C_triangulation_wi::Vertex             Vertex;
  typedef typename C_triangulation_wi::Vertex_handle      Vertex_handle;
  typedef typename C_triangulation_wi::Face_handle        Face_handle;
  typedef typename C_triangulation_wi::Edge               Edge;
  typedef typename C_triangulation_wi::Finite_faces_iterator
                                                        Finite_faces_iterator;
  typedef typename C_triangulation_wi::Face_circulator  Face_circulator;
  typedef typename CD_triangulation::Less_edge          Less_edge;
  typedef typename CD_triangulation::Edge_set           Edge_set;
  
 
  typedef std::list<Edge> List_edges;
  typedef std::list<Vertex_handle> List_vertices;  
  typedef std::list<Face_handle> List_faces;

  typedef typename Geom_traits::Point_2  Point;

  Constrained_Delaunay_triangulation_wi_2(const Geom_traits& gt=Geom_traits()) 
    : CD_triangulation(gt) { }

  Constrained_Delaunay_triangulation_wi_2(const CD_triangulation& cdt)
    : CD_triangulation(cdt) {}

  Constrained_Delaunay_triangulation_wi_2(std::list<Constraint>& lc, 
					  const Geom_traits& gt=Geom_traits())
      : CD_triangulation(gt)
  {
    typename std::list<Constraint>::iterator itc;
    itc=lc.begin();      
      do{
	insert((*itc).first, (*itc).second);
	++itc;
      } while (itc != lc.end());
      CGAL_triangulation_postcondition( is_valid() );
  }

  template<class InputIterator>
  Constrained_Delaunay_triangulation_wi_2(InputIterator first,
			      InputIterator last,
			      const Geom_traits& gt=Geom_traits() )
      : CD_triangulation(gt)
  {
    while(first != last){
          insert((*first).first, (*first).second);
	  first++;
    }
    CGAL_triangulation_postcondition( is_valid() );
  }


  // INSERTION-REMOVAL
  Vertex_handle insert(const Point & a);
  Vertex_handle special_insert_in_edge(const Point & a, Face_handle f, int i);
  void insert(const Point & a, const Point & b);
  void insert(Vertex_handle va, Vertex_handle & vb);
  void insert(Vertex_handle va, Vertex_handle vb, Face_handle & fr,
	      int & i);
  
  void flip_around(Vertex_handle va);
  void flip_around(List_vertices& new_vertices);
  void remove_2D(Vertex_handle v );

};



template < class Ctr >
inline void 
Constrained_Delaunay_triangulation_wi_2<Ctr>::
flip_around(Vertex_handle va)
  // makes the triangles incident to vertex va Delaunay using flips
{
  CD_triangulation::flip_around(va);
}

template < class Ctr >
inline void 
Constrained_Delaunay_triangulation_wi_2<Ctr>::
flip_around(List_vertices& new_vertices)
{
  typename List_vertices::iterator itv=new_vertices.begin();
  for( ; itv != new_vertices.end(); itv++) {
    flip_around(*itv);
  }
  return;
}

  
template < class Ctr >  
inline 
Constrained_Delaunay_triangulation_wi_2<Ctr>::Vertex_handle 
Constrained_Delaunay_triangulation_wi_2<Ctr>::
insert(const Point & a)
  // inserts a in the triangulation
{
  return  CD_triangulation::insert(a);
}


template < class Ctr >  
inline void 
Constrained_Delaunay_triangulation_wi_2<Ctr>::
insert(const Point & a, const Point & b)
 // inserts segment ab as a constraint and updates the 
 // constrained Delaunay triangulation
{
  Vertex_handle va= insert(a);
  Vertex_handle vb= insert(b);
  insert(va, vb);
}

template < class Ctr >  
inline void 
Constrained_Delaunay_triangulation_wi_2<Ctr>::
insert(Vertex_handle va, Vertex_handle & vb)
// inserts line segment ab as an edge in the triangulation 
{
  List_edges new_edges;
  List_vertices new_vertices;
  Face_handle fr;
  int i;
  C_triangulation_wi::insert(va,vb,fr,i,new_edges,new_vertices);
  propagating_flip(new_edges);
  flip_around(new_vertices);
}

template < class Ctr >  
inline void 
Constrained_Delaunay_triangulation_wi_2<Ctr>::
insert(Vertex_handle va, Vertex_handle vb,
       Face_handle & fr, int & i)
 // inserts line segment ab as an edge in the triangulation 
 // returns the triangle fr right to edge ab
 // edge ab=(fr,i)
{
  List_edges new_edges;
  List_vertices new_vertices;
  C_triangulation_wi::insert(va,vb,fr,i,new_edges,new_vertices);
  propagating_flip(new_edges);
  flip_around(new_vertices);
}



CGAL_END_NAMESPACE
#endif // CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_WI_2_H
