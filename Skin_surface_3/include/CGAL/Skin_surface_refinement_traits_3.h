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

#ifndef CGAL_SKIN_SURFACE_SUBDIVISION_POLICY_H
#define CGAL_SKIN_SURFACE_SUBDIVISION_POLICY_H

#include <CGAL/intersection_3_1.h>
CGAL_BEGIN_NAMESPACE

template <class Polyhedron_3, class SkinSurface_3>
class Skin_surface_subdivision_policy_base_3 {
public:
  typedef Polyhedron_3                            Polyhedron;
  typedef SkinSurface_3                           Skin_surface;
  typedef typename Polyhedron::Traits             P_traits;

  typedef typename Polyhedron::Vertex_handle      P_vertex_handle;

  typedef typename P_traits::RT          P_rt;
  typedef typename P_traits::Point_3     P_point;
  typedef typename P_traits::Segment_3   P_segment;
  typedef typename P_traits::Line_3      P_line;
  typedef typename P_traits::Vector_3    P_vector;
  typedef typename P_traits::Plane_3     P_plane;

  Skin_surface_subdivision_policy_base_3(Skin_surface const& skin)
    : ss_3(skin)
  {}
    
  virtual P_point to_surface(P_vertex_handle vh) const = 0;

  virtual P_vector normal(P_vertex_handle vh) const = 0;

protected:
  Skin_surface const &ss_3;
};

template <class Polyhedron_3, class SkinSurface_3>
class Skin_surface_subdivision_policy_default_3 
  : public Skin_surface_subdivision_policy_base_3<Polyhedron_3,SkinSurface_3>
{
public:
  typedef Polyhedron_3                            Polyhedron;
  typedef SkinSurface_3                           Skin_surface;
  typedef Skin_surface_subdivision_policy_base_3<Polyhedron_3,SkinSurface_3>
                                                  Base;
  typedef typename Polyhedron::Traits             P_traits;

  typedef typename Polyhedron::Vertex_handle      P_vertex_handle;

  typedef typename P_traits::RT          P_rt;
  typedef typename P_traits::Point_3     P_point;
  typedef typename P_traits::Segment_3   P_segment;
  typedef typename P_traits::Line_3      P_line;
  typedef typename P_traits::Vector_3    P_vector;
  typedef typename P_traits::Plane_3     P_plane;

  Skin_surface_subdivision_policy_default_3(Skin_surface const& skin)
    : Base(skin)
  {
    
  }
    
  P_point to_surface(P_vertex_handle vh) const
  {
    typename Skin_surface::Bare_point result =
      Cartesian_converter<P_traits, 
      typename Skin_surface::Geometric_traits::Kernel>()(vh->point());
    Base::ss_3.intersect_with_transversal_segment(result);
    return 
      Cartesian_converter
      <typename Skin_surface::Geometric_traits::Kernel, P_traits>()( result );
  }

  P_vector normal(P_vertex_handle vh) const
  {
    // Convert to and from the skin surface kernel
    return 
      Cartesian_converter
      <typename Skin_surface::Geometric_traits::Kernel, 
      P_traits>()( Base::ss_3.normal
		    (Cartesian_converter<P_traits, 
		     typename Skin_surface::Geometric_traits::Kernel>()
		     (vh->point())));
  }

};

template <class Polyhedron_3,
	  class SkinSurface_3>
Skin_surface_subdivision_policy_base_3<Polyhedron_3, SkinSurface_3> *
get_subdivision_policy(Polyhedron_3  &p,
		       SkinSurface_3 &skinsurface) 
{
  typedef Skin_surface_subdivision_policy_default_3<Polyhedron_3, SkinSurface_3>
    Policy;
  return new Policy(skinsurface);
}

CGAL_END_NAMESPACE

// Partial specialisation for Skin_surface_polyhedral_items_3
#include <CGAL/Skin_surface_polyhedral_items_3.h>
#include <CGAL/Skin_surface_refinement_traits_with_face_info_3.h>

CGAL_BEGIN_NAMESPACE

template <class P_Traits,
	  class SkinSurface_3>
Skin_surface_subdivision_policy_base_3<Polyhedron_3<P_Traits, 
			 Skin_surface_polyhedral_items_3<SkinSurface_3> >, SkinSurface_3> *
get_subdivision_policy(Polyhedron_3<P_Traits, 
			 Skin_surface_polyhedral_items_3<SkinSurface_3> > &p,
		       SkinSurface_3 &skinsurface) 
{
  typedef Polyhedron_3<P_Traits, 
    Skin_surface_polyhedral_items_3<SkinSurface_3> >           Polyhedron;

  typedef
    Skin_surface_subdivision_policy_with_face_info_3<Polyhedron, SkinSurface_3>
    Policy;
  
  return new Policy(skinsurface);
}

CGAL_END_NAMESPACE


#endif // CGAL_SKIN_SURFACE_SUBDIVISION_TRAITS_H
