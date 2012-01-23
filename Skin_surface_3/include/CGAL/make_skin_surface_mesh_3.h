// Copyright (c) 2005 Rijksuniversiteit Groningen (Netherlands)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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

#ifndef CGAL_MAKE_SKIN_SURFACE_MESH_3_H
#define CGAL_MAKE_SKIN_SURFACE_MESH_3_H

#include <CGAL/Skin_surface_traits_3.h>
#include <CGAL/Skin_surface_3.h>
#include <CGAL/mesh_skin_surface_3.h>
#include <CGAL/subdivide_skin_surface_mesh_3.h>

#include <CGAL/make_union_of_balls_3.h>

namespace CGAL {

template <class WP_iterator,
	  class Polyhedron_3>
void make_skin_surface_mesh_3(Polyhedron_3 &p, 
			      WP_iterator begin, WP_iterator end, 
			      double shrink_factor=.5, 
			      int nSubdivisions=0, 
			      bool grow_balls=true)
{
  if (shrink_factor == 1) {
    make_union_of_balls_mesh_3(p,begin,end,nSubdivisions);
  }


  typedef typename WP_iterator::value_type              Weighted_point;
  typedef typename Kernel_traits<Weighted_point>::Kernel K;
  
  typedef Skin_surface_traits_3<K>                      Traits;
  typedef Skin_surface_3<Traits>                        Skin_surface;
  
  Skin_surface skin_surface(begin, end, shrink_factor, grow_balls);

  mesh_skin_surface_3(skin_surface, p);

  subdivide_skin_surface_mesh_3(skin_surface, p, nSubdivisions);
}



} //namespace CGAL

#endif // CGAL_MAKE_SKIN_SURFACE_MESH_3_H
