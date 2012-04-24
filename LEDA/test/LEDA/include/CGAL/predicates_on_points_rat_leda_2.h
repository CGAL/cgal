// Copyright (c) 1999  
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
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stefan Schirra
 

#ifndef CGAL_PREDICATES_ON_POINTS_RAT_LEDA_2_H
#define CGAL_PREDICATES_ON_POINTS_RAT_LEDA_2_H

#include <CGAL/rat_leda.h>
#include <CGAL/predicates_on_points_2.h>

namespace CGAL {


template <class use_rat_leda_kernel>
inline
Orientation
orientation(const Point_2<use_rat_leda_kernel>& p,
            const Point_2<use_rat_leda_kernel>& q,
            const Point_2<use_rat_leda_kernel>& r)
{
  typedef typename use_rat_leda_kernel::Point_2_base  RPoint_2;
  return static_cast<Orientation>( CGAL_LEDA_SCOPE::orientation((const RPoint_2&)p,
                                                 (const RPoint_2&)q,
                                                 (const RPoint_2&)r) );
/*
  return static_cast<Orientation>(
      ::orientation((const use_rat_leda_kernel::Point_2_base&)p,
                    (const use_rat_leda_kernel::Point_2_base&)q,
                    (const use_rat_leda_kernel::Point_2_base&)r) );
*/
}

template <class use_rat_leda_kernel>
inline
bool
collinear(const Point_2<use_rat_leda_kernel>& p,
          const Point_2<use_rat_leda_kernel>& q,
          const Point_2<use_rat_leda_kernel>& r)
{

  typedef typename use_rat_leda_kernel::Point_2_base  RPoint_2;
  return CGAL_LEDA_SCOPE::collinear((const RPoint_2&)p,
                     (const RPoint_2&)q,
                     (const RPoint_2&)r );
/*
  return
      ::collinear((const use_rat_leda_kernel::Point_2_base&)p,
                  (const use_rat_leda_kernel::Point_2_base&)q,
                  (const use_rat_leda_kernel::Point_2_base&)r);
*/
}

} //namespace CGAL

#endif // CGAL_PREDICATES_ON_POINTS_RAT_LEDA_2_H
