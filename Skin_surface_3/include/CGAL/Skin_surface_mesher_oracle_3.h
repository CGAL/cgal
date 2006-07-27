// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Nico Kruithof

#ifndef CGAL_SKIN_SURFACE_MESHER_ORACLE_3_H
#define CGAL_SKIN_SURFACE_MESHER_ORACLE_3_H

#include <CGAL/point_generators_3.h>
//#include <CGAL/Surface_mesher/Sphere_oracle_3.h>
#include <CGAL/Surface_mesher/Oracles/Sphere_oracle_3.h>

#include <queue>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_skin_surface_3.h>

#ifdef CGAL_SURFACE_MESHER_DEBUG_CLIPPED_SEGMENT
#include <string>
#endif

#include <sstream>
#include <iostream>

// NGHK: remove later, for debugging
#include <CGAL/IO/Polyhedron_iostream.h>

CGAL_BEGIN_NAMESPACE

  template <
    class GT,
    class Surface,
    class Point_creator = Creator_uniform_3<typename GT::FT,
                                            typename GT::Point_3>
    >
  class Skin_surface_mesher_oracle_3
  {
    // private types
    typedef Skin_surface_mesher_oracle_3<GT, Surface, Point_creator> Self;

    typedef typename Surface_mesher::Sphere_oracle_3<GT, Point_creator> Sphere_oracle;
    
    typedef typename GT::Point_3 Point;

    typedef typename GT::FT FT;
    typedef typename GT::Sphere_3 Sphere_3;

  public:

    // Public types
    typedef GT Geom_traits;
    typedef typename GT::Point_3 Point_3;
    typedef typename GT::Segment_3 Segment_3;
    typedef typename GT::Ray_3 Ray_3;
    typedef typename GT::Line_3 Line_3;

    typedef Surface Surface_3;

  public:

    // Constructors
    Skin_surface_mesher_oracle_3 () 
    {
#ifdef CGAL_SURFACE_MESHER_DEBUG_CONSTRUCTORS
      std::cerr << "CONS: Skin_surface_mesher_oracle_3\n";
#endif
    }

    class Intersect_3 
    {

    public:
      Intersect_3()
      {
      }

      Object operator()(const Surface_3& surface, Segment_3 s) const
      // s is passed by value, because it is clipped below
      {
        typename GT::Construct_point_on_3 point_on =
          GT().construct_point_on_3_object();

        typename Sphere_oracle::Intersect_3 clip =
          Sphere_oracle().intersect_3_object();

        const Sphere_3& sphere = surface.bounding_sphere();

        Point_3 a = point_on(s, 0);
        Point_3 b = point_on(s, 1);

        // if both extremities are on the same side of the surface, return
        // no intersection
        if(surf_equation(surface, a) * surf_equation(surface, b) > 0)
          return Object();

        // Code for surfaces with boundaries
        // First rescale segment s = [a, b]
        if( clip.clip_segment(sphere, a, b) )
          return intersect_clipped_segment(surface,
                                           a,
                                           b,
                                           surface.squared_error_bound());
        // else
        return Object();
      } // end operator()(Surface_3, Segment_3)

      Object operator()(const Surface_3& surface, const Ray_3& r) const {
        typename Sphere_oracle::Intersect_3 clip =
          Sphere_oracle().intersect_3_object();

        const Sphere_3& sphere = surface.bounding_sphere();

        Point_3 a, b;
        if(clip.clip_ray(sphere, r, a, b))
        {
          return intersect_clipped_segment(surface,
                                           a,
                                           b,
                                           surface.squared_error_bound());
        }
        // else
        return Object();
      } // end operator()(Surface_3, Ray_3)

      Object operator()(const Surface_3& surface, const Line_3& l) const {
        typename Sphere_oracle::Intersect_3 clip =
          Sphere_oracle().intersect_3_object();

        const Sphere_3& sphere = surface.bounding_sphere();

        Point_3 a, b;
        if(clip.clip_line(sphere, l, a, b))
        {
          return intersect_clipped_segment(surface,
                                           a,
                                           b,
                                           surface.squared_error_bound());
        }
        else
          return Object();
      }; // end operator()(Surface_3, Line_3)

      // debug function
      static std::string debug_point(const Surface_3& surface,
                                     const Point& p)
      {
        std::stringstream s;
        s << p << " (distance=" 
          << CGAL::sqrt(CGAL::squared_distance(p,
                                 surface.bounding_sphere().center()))
          << ", sign=" << surf_equation(surface, p)
          << ")";
        return s.str();
      }

      static CGAL::Sign surf_equation (Surface_3 surface,
                                       const Point& p)
      {
        return CGAL::sign(surface(p));
      } // @TODO, @WARNING: we use x(), y() and z()

    private:
      // Private functions
      Object intersect_clipped_segment(const Surface_3& surface,
                                       Point p1,
                                       Point p2,
                                       const FT& squared_distance_bound) const
      {
#ifdef CGAL_SURFACE_MESHER_DEBUG_CLIPPED_SEGMENT
        std::cerr << "clipped_segment\n";
#endif
        typename GT::Compute_squared_distance_3 squared_distance = 
          GT().compute_squared_distance_3_object();
        typename GT::Construct_midpoint_3 midpoint =
          GT().construct_midpoint_3_object();

        Sign sign_at_p1 = surf_equation(surface, p1);
        Sign sign_at_p2 = surf_equation(surface, p2);

        if( sign_at_p1 == ZERO )
        {
          return make_object(p1);
        }
        if( sign_at_p2 == ZERO )
        {
          return make_object(p2);
        }

        // if both extremities are on the same side of the surface, return
        // no intersection
        if(sign_at_p1 * sign_at_p2 > 0)
          return Object();

        while(true)
        {
#ifdef CGAL_SURFACE_MESHER_DEBUG_CLIPPED_SEGMENT
          std::cerr << debug_point(surface, p1) << ", "
                    << debug_point(surface, p2) << "\n";
#endif
          Point mid = midpoint(p1, p2);
          const Sign sign_at_mid = surf_equation(surface, mid);

          if ( sign_at_mid == ZERO || 
               squared_distance(p1, p2) < squared_distance_bound )
          // If the two points are close, then we must decide
          {
#ifdef CGAL_SURFACE_MESHER_DEBUG_CLIPPED_SEGMENT
            std::cerr << "=" << debug_point(surface, mid) << "\n";
#endif
            return make_object(mid);
          }

          // Else we must go on
          if ( sign_at_p1 * sign_at_mid < 0 )
          {
            p2 = mid;
            sign_at_p2 = sign_at_mid;
          }
          else
          {
            p1 = mid;
            sign_at_p1 = sign_at_mid;
          }
        }
      } // end intersect_clipped_segment

    }; // end nested class Intersect_3

    class Construct_initial_points
    {
      Self& oracle;
    public:
      Construct_initial_points(Self& oracle) : oracle(oracle)
      {
      }
      
      // Random points
      template <typename OutputIteratorPoints>
      OutputIteratorPoints operator() (const Surface_3& surface, 
                                       OutputIteratorPoints out, 
                                       int n = 20) // WARNING: why 20?
      {
	typedef CGAL::Polyhedron_3<GT> Polyhedron;
	typedef typename Polyhedron::Vertex_iterator V_iterator;
	typedef typename Polyhedron::Halfedge_around_vertex_circulator 
	                                               HAV_circulator;

	Polyhedron coarse_mesh;
	mesh_skin_surface_3(surface, coarse_mesh);
	
	for (V_iterator vit = coarse_mesh.vertices_begin();
	     vit != coarse_mesh.vertices_end(); ) {
	  HAV_circulator hav_start, hav_it;
	  hav_start = hav_it = vit->vertex_begin();
	  bool ok = true;
	  FT dist = surface.get_density(vit->point());
	  do {
	    ok = (squared_distance(vit->point(),
				   hav_it->vertex()->point()) < dist);
	  } while (ok && (hav_start != (++hav_it)));
	  
	  vit++;
	  if (!ok) {
	    std::cout << "remove Vertex" << std::endl;
	    V_iterator tmp = vit; tmp--;
	    CGAL_assertion(squared_distance(tmp->point(),
					    hav_it->vertex()->point()) < dist);
	    if (tmp->degree() == 3) {
	      CGAL_assertion(hav_it->opposite()->vertex() == tmp);
	      coarse_mesh.erase_center_vertex(hav_it->opposite());
	      continue;
	    }
	    if (hav_it->vertex()->degree() == 3) {
	      coarse_mesh.erase_center_vertex(hav_it);
	      continue;
	    }
	    coarse_mesh.join_vertex(hav_it);
	  }
	}

	std::copy(coarse_mesh.points_begin(), coarse_mesh.points_end(), out);
	{
	  std::ofstream os("delaunay_coarsed_mesh.off");
	  os << coarse_mesh;
	}
        return out;
      }
    }; // end nested class Construct_initial_points

    Construct_initial_points construct_initial_points_object()
    {
      return Construct_initial_points(*this);
    }

    Intersect_3 intersect_3_object()
    {
      return Intersect_3();
    }

    bool is_in_volume(const Surface_3& surface, const Point& p)
    {
      return Intersect_3::surf_equation(surface, p)<0.;
    }

  };  // end Skin_surface_mesher_oracle_3

template <typename FT>
FT approximate_sqrt(const FT x) {
  return FT (CGAL_NTS sqrt(CGAL_NTS to_double(x)));
}

CGAL_END_NAMESPACE


#endif  // CGAL_SKIN_SURFACE_MESHER_ORACLE_3_H
