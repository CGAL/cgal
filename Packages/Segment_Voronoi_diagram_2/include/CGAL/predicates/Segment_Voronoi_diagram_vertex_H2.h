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




#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_VORONOI_VERTEX_H2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_VORONOI_VERTEX_H2_H


#include <CGAL/Number_type_traits.h>
//#include <CGAL/predicates/Svd_Voronoi_vertex_ring_H2.h>
#include <CGAL/predicates/Svd_Voronoi_vertex_sqrt_field_H2.h>


CGAL_BEGIN_NAMESPACE



template<class M>
struct Svd_which_base_H2
{
  template<class K>
  struct Which_base {
    // this should change once I have written the code for the ring
    //    operations for a homogeneous kernel
    //    typedef Svd_voronoi_vertex_ring_H2<K>  Base;
    typedef Svd_voronoi_vertex_sqrt_field_H2<K>  Base;
  };
};

template <>
struct Svd_which_base_H2<Sqrt_field_tag>
{
  template <class K>
  struct Which_base {
    typedef Svd_voronoi_vertex_sqrt_field_H2<K>  Base;
  };
};


//----------------------------------------------------------------------

template<class K, class M>
class Svd_voronoi_vertex_H2
  : public Svd_which_base_H2<M>::template Which_base<K>::Base
{
private:
  typedef typename
  Svd_which_base_H2<M>::template Which_base<K>::Base  Base;

protected:
  typedef typename Base::Site_2   Site_2;
public:
  Svd_voronoi_vertex_H2(const Site_2& p, const Site_2& q,
			const Site_2& r)
    : Base(p, q, r) {}
};



CGAL_END_NAMESPACE



#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_VORONOI_VERTEX_H2_H
