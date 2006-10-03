// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
#ifndef CGAL_TYPEDEFS_H
#define CGAL_TYPEDEFS_H

#include <CGAL/Cartesian.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/envelope_3.h>
#include <CGAL/Env_triangle_traits_3.h>
#include <CGAL/Bbox_2.h> 

#ifdef CGAL_USE_GMP

  #include <CGAL/Gmpq.h>

  typedef CGAL::Gmpq                                    Base_nt;

#else

  #include <CGAL/MP_Float.h>
  #include <CGAL/Quotient.h>

  typedef CGAL::Quotient<CGAL::MP_Float>                Base_nt;

#endif

typedef CGAL::Lazy_exact_nt<Base_nt>                  Coord_type;
typedef CGAL::Cartesian<Coord_type>		                Kernel;

typedef Kernel::Point_2				                        Point_2;
typedef Kernel::Point_3				                        Point_3;
typedef Kernel::Triangle_2				                    Triangle_2;
typedef Kernel::Triangle_3				                    Triangle_3;


typedef CGAL::Env_triangle_traits_3<Kernel>           Traits;
typedef Traits::Surface_3                             Surface_3;
typedef Traits::Xy_monotone_surface_3                 Xy_monotone_surface_3;
typedef CGAL::Envelope_diagram_2<Traits>              Envelope_diagram_2;
typedef Envelope_diagram_2::Edge_const_iterator       Edge_const_iterator;
typedef Envelope_diagram_2::Vertex_const_iterator     Vertex_const_iterator;

CGAL::Bbox_2 bbox_2d(const Triangle_3& tri)
{
  const Point_3& p0 = tri.vertex(0);
  const Point_3& p1 = tri.vertex(1);
  const Point_3& p2 = tri.vertex(2);
  Triangle_2 proj_tri(Point_2(p0.x(), p0.y()),
                      Point_2(p1.x(), p1.y()),
                      Point_2(p2.x(), p2.y()));
  return proj_tri.bbox();
}
#endif
