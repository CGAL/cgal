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
// file          : include/CGAL/Constrained_Delaunay_triangulation_wi_2.h
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
template < class  Ctwi_ >
class Constrained_Delaunay_triangulation_wi_base_2;

template <class Gt, class Tds>
class Constrained_Delaunay_triangulation_wi_2
  : public Constrained_Delaunay_triangulation_wi_base_2<
           Constrained_triangulation_wi_2<Gt, Tds> >
{
public:
  typedef Constrained_triangulation_wi_2<Gt, Tds>               Ctwi;
  typedef Constrained_Delaunay_triangulation_wi_base_2<Ctwi>    Base;
  typedef Constrained_Delaunay_triangulation_wi_2<Gt, Tds>      CDtwi;
  typedef typename Base::Geom_traits   Geom_traits;
  typedef typename Base::Constraint    Constraint;

  Constrained_Delaunay_triangulation_wi_2(const Geom_traits& gt=Geom_traits()) 
    : Base(gt) { }

  Constrained_Delaunay_triangulation_wi_2(const CDtwi& cdt)
    : Base(cdt) {}

  Constrained_Delaunay_triangulation_wi_2(std::list<Constraint>& lc, 
				       const Geom_traits& gt=Geom_traits())
    : Base(lc,gt) {}

  template<class InputIterator>
  Constrained_Delaunay_triangulation_wi_2(InputIterator first,
				       InputIterator last,
				       const Geom_traits& gt=Geom_traits() )
    : Base(first,last,gt) {}

}; 

template < class  Ctwi_ >
class Constrained_Delaunay_triangulation_wi_base_2
  : public Constrained_Delaunay_triangulation_base_2<Ctwi_>
{
public:
  typedef Ctwi_                                               Ctwi;
  typedef Constrained_Delaunay_triangulation_base_2<Ctwi>     CDt_base;
  typedef Constrained_Delaunay_triangulation_wi_base_2<Ctwi>  CDtwi_base;
  typedef typename Ctwi::Geom_traits        Geom_traits;
  typedef typename Ctwi::Constraint         Constraint;
  typedef typename Ctwi::Vertex             Vertex;
  typedef typename Ctwi::Vertex_handle      Vertex_handle;
  typedef typename Ctwi::Face_handle        Face_handle;
  typedef typename Ctwi::Edge               Edge;
  typedef typename Ctwi::Finite_faces_iterator      Finite_faces_iterator;
  typedef typename Ctwi::Face_circulator            Face_circulator;
  typedef typename CDt_base::Less_edge              Less_edge;
  typedef typename CDt_base::Edge_set               Edge_set;
 
  typedef std::list<Edge> List_edges;
  typedef std::list<Vertex_handle> List_vertices;  
  typedef std::list<Face_handle> List_faces;

  typedef typename Geom_traits::Point_2  Point;

  Constrained_Delaunay_triangulation_wi_base_2(const Geom_traits& 
					       gt=Geom_traits()) 
    : CDt_base(gt) { }

  Constrained_Delaunay_triangulation_wi_base_2(const CDtwi_base& cdt)
    : CDt_base(cdt) {}

  Constrained_Delaunay_triangulation_wi_base_2(std::list<Constraint>& lc, 
					  const Geom_traits& gt=Geom_traits())
      : CDt_base(gt)
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
  Constrained_Delaunay_triangulation_wi_base_2(InputIterator first,
					       InputIterator last,
			      const Geom_traits& gt=Geom_traits() )
      : CDt_base(gt)
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
Constrained_Delaunay_triangulation_wi_base_2<Ctr>::
flip_around(Vertex_handle va)
  // makes the triangles incident to vertex va Delaunay using flips
{
  CDt_base::flip_around(va);
}

template < class Ctr >
inline void 
Constrained_Delaunay_triangulation_wi_base_2<Ctr>::
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
Constrained_Delaunay_triangulation_wi_base_2<Ctr>::Vertex_handle 
Constrained_Delaunay_triangulation_wi_base_2<Ctr>::
insert(const Point & a)
  // inserts a in the triangulation
{
  return  CDt_base::insert(a);
}


template < class Ctr >  
inline void 
Constrained_Delaunay_triangulation_wi_base_2<Ctr>::
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
Constrained_Delaunay_triangulation_wi_base_2<Ctr>::
insert(Vertex_handle va, Vertex_handle & vb)
// inserts line segment ab as an edge in the triangulation 
{
  List_edges new_edges;
  List_vertices new_vertices;
  Face_handle fr;
  int i;
  Ctwi::insert(va,vb,fr,i,new_edges,new_vertices);
  propagating_flip(new_edges);
  flip_around(new_vertices);
}

template < class Ctr >  
inline void 
Constrained_Delaunay_triangulation_wi_base_2<Ctr>::
insert(Vertex_handle va, Vertex_handle vb,
       Face_handle & fr, int & i)
 // inserts line segment ab as an edge in the triangulation 
 // returns the triangle fr right to edge ab
 // edge ab=(fr,i)
{
  List_edges new_edges;
  List_vertices new_vertices;
  Ctwi::insert(va,vb,fr,i,new_edges,new_vertices);
  propagating_flip(new_edges);
  flip_around(new_vertices);
}



CGAL_END_NAMESPACE
#endif // CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_WI_2_H
