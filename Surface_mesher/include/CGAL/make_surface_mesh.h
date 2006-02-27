// Copyright (c) 2006  INRIA Sophia-Antipolis (France).
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
// $Source: 
// $Revision: 1.1 $ $Date: 2005/12/12 16:20:58 $
// $Name:  $
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_MAKE_SURFACE_MESH_H
#define CGAL_MAKE_SURFACE_MESH_H

#include <CGAL/Surface_mesher/Surface_mesher.h
#include <CGAL/Surface_mesh_traits_3_generator.h>

namespace CGAL {

template <typename C2T3,
	  typename Surface,
	  typename Criteria,
	  typename Tag>
void make_surface_mesh(C2T3& c2t3,
                       Surface surface,
                       Criteria criteria,
                       Tag tag = Non_manifold_tag() ) 
{
  typedef typename CGAL::Surface_mesh_traits_3_generator<Surface>::type Traits;

  make_surface_mesh(c2t3, surface, criteria, Traits(), tag);  
}

template <typename C2T3,
	  typename SurfaceMeshTraits_3,
	  typename Criteria>
void make_surface_mesh(C2T3& c2t3,
                       typename SurfaceMeshTraits_3::Surface_3 surface,
                       Criteria criteria,
		       SurfaceMeshTraits_3 traits,
                       Non_manifold_tag)
{
  using namespace CGAL::Surface_mesher;
  
  typedef Surface_mesher<
    typename SurfaceMeshTraits_3::Surface_3,
    Criteria,
    C2T3> Mesher;
  
  Mesher mesher(c2t3, surface, criteria, traits);
  // TODO initial, then refine()
}

} // end namespace CGAL

#endif CGAL_MAKE_SURFACE_MESH_H
