// Copyright (c) 2003  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Philippe Guigue

#ifndef CGAL_TRIANGLE_3_PLANE_3_DO_INTERSECT_H
#define CGAL_TRIANGLE_3_PLANE_3_DO_INTERSECT_H


CGAL_BEGIN_NAMESPACE

namespace CGALi {

template <class K>
bool do_intersect(const typename CGAL_WRAP(K)::Triangle_3 &t, 
		  const typename CGAL_WRAP(K)::Plane_3   &h,
		  const K & k)
{
  
  CGAL_kernel_precondition( ! k.is_degenerate_3_object()(t)) ;
  CGAL_kernel_precondition( ! k.is_degenerate_3_object()(h)) ;

   
  typename K::Construct_vertex_3 vertex_on =
    k.construct_vertex_3_object();
  
  typename K::Oriented_side_3 oriented_side =
    k.oriented_side_3_object();

  

  switch ( oriented_side(h,vertex_on(t,0)) ) {
  case ON_POSITIVE_SIDE:
    return oriented_side(h,vertex_on(t,1)) != ON_POSITIVE_SIDE
      || oriented_side(h,vertex_on(t,2)) != ON_POSITIVE_SIDE;
  case ON_NEGATIVE_SIDE:
    return oriented_side(h,vertex_on(t,1)) != ON_NEGATIVE_SIDE
      || oriented_side(h,vertex_on(t,2)) != ON_NEGATIVE_SIDE ;
  case ON_ORIENTED_BOUNDARY:
    return true;
  default:// should not happen.
    CGAL_kernel_assertion(false);
    return false;
  }
}


template <class K>
inline
bool do_intersect(const typename CGAL_WRAP(K)::Plane_3   &h,
		  const typename CGAL_WRAP(K)::Triangle_3 &t, 
		  const K & k)
{
  return do_intersect(t, h, k);
}


} // namespace CGALi


template <class K>
inline bool do_intersect(const Triangle_3<K> &t, 
			 const Plane_3<K>    &h)
{
  return typename K::Do_intersect_3()(t,h);
}
  
template <class K>
inline bool do_intersect(const Plane_3<K>    &h, 
			 const Triangle_3<K> &t)
{
  return typename K::Do_intersect_3()(t,h);
}

/*
template <class K>
inline bool do_intersect(const Plane_3<K>    &h, 
			 const Triangle_3<K> &t,
			 const K & k)
{
  return do_intersect(t,h,k);
}
*/


CGAL_END_NAMESPACE


#endif //CGAL_TRIANGLE_3_PLANE_3_DO_INTERSECT_H




