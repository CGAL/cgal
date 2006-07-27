// Copyright (c) 2005 Rijksuniversiteit Groningen (Netherlands)
// All rights reserved.
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
// $URL$
// $Id$
// 
//
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef CGAL_MIXED_COMPLEX_TRAITS_3_H
#define CGAL_MIXED_COMPLEX_TRAITS_3_H

#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
// Contains the weighted converter:
#include <CGAL/Regular_triangulation_filtered_traits_3.h>

CGAL_BEGIN_NAMESPACE 

/** Input: a list of n weighted points p_1...p_n and a query point x.
    There is a plane separating the mixed cell defined by p_1...p_n-1
    and the mixed cell defined by p_1...p_n. The predicate tests
    whether x lies on the same side of this plane as the mixed cell
    defined by p_1...p_n-1 (NEGATIVE), on the plane (ZERO) or on the
    opposite side (POSITIVE).
 **/
template <class K>
class Side_of_mixed_cell {
public:
  typedef typename K::FT             FT;
  typedef typename K::Bare_point     Bare_point;
  typedef typename K::Weighted_point Weighted_point;

  Side_of_mixed_cell(const FT &shrink) : s(shrink) {}
  
  typedef CGAL::Arity_tag< 5 > Arity;
  typedef CGAL::Sign           result_type;
  
  Test_add() : _t(1) {}
  Test_add(FT t) : _t(t) {}
  
  result_type operator()(const Weighted_point &p1,
			 const Weighted_point &p2,
			 const Bare_point &x) const {
    return side_of_mixed_cellC3(p1.x(),p1.y(),p1.z(),p1.weight(),
				p2.x(),p2.y(),p2.z(),p2.weight(),
				x.x(),x.y(),x.z(),
				s);
  }
  result_type operator()(const Weighted_point &p1,
			 const Weighted_point &p2,
			 const Weighted_point &p3,
			 const Bare_point &x) const {
    return side_of_mixed_cellC3(p1.x(),p1.y(),p1.z(),p1.weight(),
				p2.x(),p2.y(),p2.z(),p2.weight(),
				p3.x(),p3.y(),p3.z(),p3.weight(),
				x.x(),x.y(),x.z(),
				s);
  }
  result_type operator()(const Weighted_point &p1,
			 const Weighted_point &p2,
			 const Weighted_point &p3,
			 const Weighted_point &p4,
			 const Bare_point &x) const {
    return side_of_mixed_cellC3(p1.x(),p1.y(),p1.z(),p1.weight(),
				p2.x(),p2.y(),p2.z(),p2.weight(),
				p3.x(),p3.y(),p3.z(),p3.weight(),
				p4.x(),p4.y(),p4.z(),p4.weight(),
				x.x(),x.y(),x.z(),
				s);
  }
  
private:
  FT s;
}

template <class K>
class Mixed_complex_traits_3_base {
  
};


CGAL_END_NAMESPACE

#endif // CGAL_MIXED_COMPLEX_TRAITS_3_H
