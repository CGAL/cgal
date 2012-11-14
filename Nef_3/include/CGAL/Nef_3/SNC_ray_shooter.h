// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
#ifndef CGAL_SNC_RAY_SHOOTER_H
#define CGAL_SNC_RAY_SHOOTER_H

#include <CGAL/basic.h>
#include <CGAL/functional.h> 
#include <CGAL/function_objects.h> 
#include <CGAL/Nef_3/Pluecker_line_3.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_SM_overlayer.h>
#include <CGAL/Nef_S2/SM_point_locator.h>
#include <CGAL/Nef_3/SNC_FM_decorator.h>
#include <CGAL/Nef_3/SNC_intersection.h>

#ifdef SM_VISUALIZOR
#include <CGAL/Nef_3/SNC_SM_visualizor.h>
#endif // SM_VISUALIZOR
#include <map>
#include <list>
#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 37
#include <CGAL/Nef_2/debug.h>

#ifndef CGAL_I_DO_WANT_TO_USE_GENINFO
#include <boost/any.hpp>
#endif

namespace CGAL {

// ----------------------------------------------------------------------------
// SNC_ray_shooting
// ----------------------------------------------------------------------------

/*{\Manpage{SNC_ray_shooting}{SNC}{ray shoot functionality}{O}}*/

template <typename SNC_decorator>
class SNC_ray_shooter : public SNC_decorator
{ 

protected:
  typedef typename SNC_decorator::SNC_structure   SNC_structure;
  typedef SNC_ray_shooter<SNC_decorator>          Self;
  typedef SNC_decorator                           Base;

public:
  typedef typename SNC_decorator::Decorator_traits Decorator_traits;
  typedef typename Decorator_traits::SM_decorator SM_decorator;
  typedef SM_point_locator<SM_decorator>           SM_point_locator;
  typedef SNC_intersection<SNC_structure>          SNC_intersection;

  typedef typename Decorator_traits::Vertex_handle Vertex_handle;
  typedef typename Decorator_traits::Halfedge_handle Halfedge_handle;
  typedef typename Decorator_traits::Halffacet_handle Halffacet_handle;
  typedef typename Decorator_traits::Volume_handle Volume_handle;

  typedef typename Decorator_traits::SVertex_handle SVertex_handle;
  typedef typename Decorator_traits::SHalfedge_handle SHalfedge_handle;
  typedef typename Decorator_traits::SFace_handle SFace_handle;
  typedef typename Decorator_traits::SHalfloop_handle SHalfloop_handle;

  typedef typename SNC_structure::Object_handle Object_handle;

  typedef typename SNC_structure::Kernel Kernel;
  typedef typename SNC_structure::Point_3 Point_3;
  typedef typename SNC_structure::Vector_3 Vector_3;
  typedef typename SNC_structure::Segment_3 Segment_3;
  typedef typename SNC_structure::Ray_3 Ray_3;
  typedef typename SNC_structure::Line_3 Line_3;
  typedef typename SNC_structure::Plane_3 Plane_3;

  typedef typename SNC_structure::Mark Mark;

  #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
  typedef void* GenPtr;
  #else
  typedef boost::any GenPtr;
  #endif

  SNC_ray_shooter() {}
  void initialize(SNC_structure* W) { *this = SNC_ray_shooter(*W);}

  SNC_ray_shooter(SNC_structure& W) : Base(W) {}
  /*{\Mcreate makes |\Mvar| a ray shooter on |W|.}*/

 private:
  Volume_handle determine_volume(const Ray_3& ray) const {
    CGAL_precondition( !ray.is_degenerate());
    Object_handle o = shoot(ray);
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f, f_below;
    if( CGAL::assign(v, o)) {
      CGAL_NEF_TRACEN("facet below from from vertex...");
      f_below = get_visible_facet(v, ray);
      if(f_below != Halffacet_handle())
	return f_below->incident_volume();
      SM_decorator SD(&*v);
      CGAL_assertion( SD.number_of_sfaces() == 1);
      return SD.sfaces_begin()->volume();
    }
    else if( CGAL::assign(e, o)) {
      CGAL_NEF_TRACEN("facet below from from edge...");
      f_below = get_visible_facet(e, ray);
      if(f_below != Halffacet_handle())
	return f_below->incident_volume();
      SM_decorator SD(&*e->source());
      CGAL_assertion(SD.is_isolated(e));
      return e->incident_sface()->volume();
    }
    else if( CGAL::assign(f, o)) {
      CGAL_NEF_TRACEN("facet below from from facet...");
      f_below = get_visible_facet(f, ray);
      CGAL_assertion( f_below != Halffacet_handle());
      return f_below->incident_volume();
    }
    
    return Base(*this).volumes_begin();
  }

 public:
  Object_handle shoot(const Ray_3& ray) const
     /*{\Mop returns the nearest object hit by a ray |ray|. }*/ {
    CGAL_precondition( !ray.is_degenerate());
    bool hit = false;
    Point_3 end_of_seg;
    SNC_intersection is(*this->sncp());

    CGAL_NEF_TRACEN( "Shooting ray " << ray);
    Object_handle o;
    Vertex_handle v;
    CGAL_forall_vertices( v, *this->sncp()) {
      if ( ray.source() != v->point() && ray.has_on(v->point())) {
        if(hit && !Segment_3(ray.source(), end_of_seg).has_on(v->point()))
          continue;
        CGAL_NEF_TRACEN("ray hit vertex case "<<v->point());
        end_of_seg = v->point();
        hit = true;
        o = Object_handle(v);
      }
    }

    Halfedge_handle e;
    CGAL_forall_edges( e, *this->sncp()) {
      Point_3 q;
      if( is.does_intersect_internally( ray, segment(e), q)) {
        if (!hit || 
	    has_smaller_distance_to_point(ray.source(),q, end_of_seg)) {
          CGAL_NEF_TRACEN("ray hit edge case " << segment(e) << " in " << q);
          end_of_seg = q;
          hit = true;
          o = Object_handle(e);
        }
      }
    }

    Halffacet_handle f;
    CGAL_forall_halffacets( f, *this->sncp()) {
      Point_3 q;
      if( is.does_intersect_internally( ray, f, q) ) {
        if(!hit || 
	   has_smaller_distance_to_point(ray.source(), q, end_of_seg)) {
        CGAL_NEF_TRACEN("ray hit facet "<< f->plane()<<" on "<<q);
        end_of_seg = q;
        hit = true;
        o = Object_handle(f);
        }
      }
    }
    return o;
  }

  Object_handle locate( const Point_3& p) const
    /*{\Mop returns the lowest dimension object on an SNC structure
      which contais |p| in its interior. }*/ {

    SNC_intersection is(*this->sncp());

    CGAL_NEF_TRACEN( "Point locator for " << p);
    Vertex_handle v;
    CGAL_forall_vertices( v, *this->sncp()) {
      CGAL_NEF_TRACEN("test vertex " << v->point());
      if ( p == v->point()) {
	CGAL_NEF_TRACEN("on vertex.");
	return Object_handle(v);
      }
    }

    Halfedge_handle e;
    CGAL_forall_edges( e, *this->sncp()) {
      if ( is.does_contain_internally( segment(e), p) ) {
	CGAL_NEF_TRACEN("on edge.");
	return Object_handle(e);
      }
    }
    Halffacet_handle f;
    CGAL_forall_halffacets( f, *this->sncp()) {
      if ( is.does_contain_internally( f, p) ) {
	CGAL_NEF_TRACEN("on facet.");
	return Object_handle(f);
      }
    }

    CGAL_warning("altered code in SNC_ray_shooter");
    Ray_3 r( p, Vector_3( 0, 0, 1));
    return Object_handle(determine_volume(r));
  }   

}; // SNC_ray_shooter

} //namespace CGAL

#endif //CGAL_SNC_RAY_SHOOTER_H
