// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France) and
// Notre Dame University (U.S.A.).  All rights reserved.
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
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>




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
