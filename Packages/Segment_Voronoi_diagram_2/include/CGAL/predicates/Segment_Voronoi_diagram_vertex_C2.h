// ======================================================================
//
// Copyright (c) 2003 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/predicates/Segment_Voronoi_diagram_vertex_C2.h
// package       : Segment_Voronoi_diagram_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   :
//
// ======================================================================



#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_VORONOI_VERTEX_C2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_VORONOI_VERTEX_C2_H


#include <CGAL/Number_type_traits.h>
#include <CGAL/predicates/Svd_Voronoi_vertex_ring_C2.h>
#include <CGAL/predicates/Svd_Voronoi_vertex_sqrt_field_C2.h>


CGAL_BEGIN_NAMESPACE



template<class M>
struct Svd_which_base_C2
{
  template<class K>
  struct Which_base {
    typedef Svd_voronoi_vertex_ring_C2<K>  Base;
  };
};

template <>
struct Svd_which_base_C2<Sqrt_field_tag>
{
  template <class K>
  struct Which_base {
    typedef Svd_voronoi_vertex_sqrt_field_C2<K>  Base;
  };
};


//----------------------------------------------------------------------

template<class K, class M>
class Svd_voronoi_vertex_2
  : public Svd_which_base_C2<M>::template Which_base<K>::Base
{
private:
  typedef typename
  Svd_which_base_C2<M>::template Which_base<K>::Base  Base;

  typedef typename Base::Site_2   Site_2;
public:
  Svd_voronoi_vertex_2(const Site_2& p, const Site_2& q,
		       const Site_2& r)
    : Base(p, q, r) {}
};



CGAL_END_NAMESPACE



#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_VORONOI_VERTEX_C2_H
