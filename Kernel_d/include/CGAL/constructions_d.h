// Copyright (c) 2000,2001  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Michael Seel

#ifndef CGAL_CONSTRUCTIONS_D_H
#define CGAL_CONSTRUCTIONS_D_H

#include <CGAL/Kernel_d/Point_d.h>

namespace CGAL {

/*{\Moptions outfile=constructions_d.man}*/

/*{\Mtext \setopdims{4cm}{2cm}\computewidths
\headerline{Constructions}}*/

template <class R>
Point_d<R> lift_to_paraboloid(const Point_d<R>& p)
/*{\Mfunc returns $p = (x_0,\ldots,x_{d-1})$ lifted to the paraboloid 
of revolution . }*/
{ typename R::Lift_to_paraboloid_d lift; 
  return lift(p); }

template <class R>
Point_d<R> project_along_d_axis(const Point_d<R>& p)
/*{\Mfunc returns |p| projected along the $d$-axis onto the hyperspace
spanned by the first $d-1$ standard base vectors.}*/
{ typename R::Project_along_d_axis_d proj; 
  return proj(p); }

template <class R>
Point_d<R> midpoint(const Point_d<R>& p, const Point_d<R>& q)
/*{\Mfunc returns the midpoint of |p| and |q|. }*/
{ typename R::Midpoint_d mp; return mp(p,q); }

template <class R, class Forward_iterator>
Point_d<R> center_of_sphere(Forward_iterator start, Forward_iterator end)
{ typename R::Center_of_sphere_d center;
  return center(start,end); }

template <class R> typename R::FT
squared_distance(const Point_d<R>& p, const Point_d<R>& q)
/*{\Mfunc returns the squared distance between |p| and $q$. }*/
{ typename R::Squared_distance_d dist; return dist(p,q); }

template <class ForwardIterator, class OutputIterator>
OutputIterator linear_base(
  ForwardIterator first, ForwardIterator last, OutputIterator result)
/*{\Mfunc computes a basis of the linear space spanned by the vectors 
in |set [first,last)| and returns it via an iterator range starting
in |result|. The returned iterator marks the end of the output.
\precond value type of |ForwardIterator| and |OutputIterator| is 
|Vector_d<R>|.}*/
{ typedef typename std::iterator_traits<ForwardIterator>::
    value_type value_type;
  typedef typename value_type::R R;
  typename R::Linear_base_d base; 
  return base(first,last,result); 
}


} //namespace CGAL
#endif // CGAL_CONSTRUCTIONS_D_H
