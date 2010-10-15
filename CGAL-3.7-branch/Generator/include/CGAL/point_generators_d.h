// Copyright (c) 1997  Utrecht University (The Netherlands),
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
// $URL$
// $Id$
// 
//
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>

#ifndef CGAL_POINT_GENERATORS_D_H
#define CGAL_POINT_GENERATORS_D_H 1

#include <CGAL/generators.h>
#include <CGAL/iterator.h>
#include <iterator>

namespace CGAL {
template < class P>
class Random_points_in_iso_box_d : public Generator_base<P> 
{
  void generate_point();
  
  typedef Counting_iterator<Random_double_in_interval> N_Random_double_iterator;
  typedef Creator_uniform_d<N_Random_double_iterator,P> Creator;
  int d;
  Random_double_in_interval rdii;

  typedef Random_points_in_iso_box_d<P> This;
 public:
  Random_points_in_iso_box_d()
    {}

  Random_points_in_iso_box_d(int dim, double a = 1, Random& rnd = default_random)
    : Generator_base<P>(a), d(dim), rdii(a, rnd)
    { 
      generate_point(); 
    }
  
  This& operator++() {
    generate_point();
    return *this;
  }
  This  operator++(int) {
    This tmp = *this;
    ++(*this);
    return tmp;
  }
};

template < class P>
void
Random_points_in_iso_box_d<P>::
generate_point() 
{
  Creator creator(d);
  rdii++;
  this->d_item =
    creator( N_Random_double_iterator(rdii), N_Random_double_iterator(d));
}

} //namespace CGAL
#endif // CGAL_POINT_GENERATORS_D_H //
// EOF //
