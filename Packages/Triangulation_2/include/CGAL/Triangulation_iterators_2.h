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
// file          : Triangulation/include/CGAL/Triangulation_iterators_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_ITERATORS_2_H
#define CGAL_TRIANGULATION_ITERATORS_2_H

#include <utility>
#include <iterator>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_utils_2.h>
#include <CGAL/Triangulation_ds_iterators_2.h>

CGAL_BEGIN_NAMESPACE

template < class Triangulation>
class Triangulation_finite_faces_iterator_2
 : public Filter_iterator< typename Triangulation::All_faces_iterator, 
                           typename Triangulation::Infinite_tester >
{
public:
  typedef typename Triangulation::All_faces_iterator        All_faces_iterator;
  typedef typename Triangulation::Infinite_tester           Infinite_tester;
  typedef Filter_iterator<All_faces_iterator,Infinite_tester>  Base;
  typedef Triangulation_finite_faces_iterator_2<Triangulation> Self;
  typedef typename Triangulation::Face_handle               Face_handle;
  
  Triangulation_finite_faces_iterator_2()    : Base() {}

  Triangulation_finite_faces_iterator_2(const Triangulation* tr)
    : Base( filter_iterator(tr->all_faces_begin(), 
			    tr->all_faces_end(),
			    tr->infinite_tester())) {}

  Triangulation_finite_faces_iterator_2(const Triangulation* tr, int)
    : Base( filter_iterator(tr->all_faces_begin(), 
			    tr->all_faces_end(),
			    tr->infinite_tester(),
			    tr->all_faces_end())) {}

  operator Face_handle() const {return (*this)->handle();}
  Self&  operator++() { Base::operator++(); return *this;}
  Self&  operator--() { Base::operator--(); return *this; }
  Self   operator++(int) { Self tmp(*this); ++(*this); return tmp; }
  Self   operator--(int) { Self tmp(*this); --(*this); return tmp; }
};


template < class Triangulation>
class Triangulation_finite_vertices_iterator_2
 : public Filter_iterator< typename Triangulation::All_vertices_iterator, 
                           typename Triangulation::Infinite_tester >
{
public:
  typedef typename Triangulation::All_vertices_iterator All_vertices_iterator;
  typedef typename Triangulation::Infinite_tester       Infinite_tester;
  typedef Filter_iterator<All_vertices_iterator,Infinite_tester>  Base;
  typedef Triangulation_finite_vertices_iterator_2<Triangulation> Self;
  typedef typename Triangulation::Vertex_handle         Vertex_handle;

  Triangulation_finite_vertices_iterator_2()    :  Base() {}
           
  Triangulation_finite_vertices_iterator_2(const Triangulation *tr)
    : Base( filter_iterator(tr->all_vertices_begin(), 
			    tr->all_vertices_end(),
			    tr->infinite_tester())) {}
    
  Triangulation_finite_vertices_iterator_2(const Triangulation *tr, int )
    : Base( filter_iterator(tr->all_vertices_begin(), 
			    tr->all_vertices_end(),
			    tr->infinite_tester(),
			    tr->all_vertices_end())) {}
			    
  operator Vertex_handle() const {return (*this)->handle();}  
  Self&  operator++() { Base::operator++(); return *this;}
  Self&  operator--() { Base::operator--(); return *this;}
  Self   operator++(int) { Self tmp(*this); ++(*this); return tmp; }
  Self   operator--(int) { Self tmp(*this); --(*this); return tmp; }

};  



CGAL_END_NAMESPACE


#endif //CGAL_TRIANGULATION_ITERATORS_2_H
