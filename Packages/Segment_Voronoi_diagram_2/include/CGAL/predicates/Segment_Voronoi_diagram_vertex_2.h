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
// file          : include/CGAL/predicates/Segment_Voronoi_diagram_vertex_2.h
// package       : Segment_Voronoi_diagram_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   :
//
// ======================================================================



#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_VORONOI_VERTEX_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_VORONOI_VERTEX_2_H


#include <CGAL/predicates/Segment_Voronoi_diagram_vertex_C2.h>
#include <CGAL/predicates/Segment_Voronoi_diagram_vertex_H2.h>


CGAL_BEGIN_NAMESPACE

template<class Rep> struct Svd_which_rep_2;

template<>
struct Svd_which_rep_2<Cartesian_tag>
{
  template<class K, class M>
  struct Which_rep {
    typedef Svd_voronoi_vertex_C2<K,M>  Rep;
  };
};

template<>
struct Svd_which_rep_2<Homogeneous_tag>
{
  template <class K, class M>
  struct Which_base {
    typedef Svd_voronoi_vertex_H2<K,M>  Rep;
  };
};


//----------------------------------------------------------------------

template<class K, class M>
class Svd_voronoi_vertex_2
  : public Svd_which_rep_2<typename K::Rep_tag>::template Which_rep<K,M>::Rep
{
private:
  typedef typename K::Rep_tag  Rep_tag;
  typedef typename
  Svd_which_rep_2<Rep_tag>::template Which_rep<K,M>::Rep  Base;

  typedef typename Base::Site_2   Site_2;
public:
  Svd_voronoi_vertex_2(const Site_2& p, const Site_2& q,
		       const Site_2& r)
    : Base(p, q, r) {}
};



CGAL_END_NAMESPACE



#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_VORONOI_VERTEX_2_H
