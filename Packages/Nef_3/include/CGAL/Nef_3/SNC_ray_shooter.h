// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
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
// $Source$
// $Revision$ $Date$
// $Name$
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
#include <CGAL/Circulator_project.h>
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
#undef _DEBUG
#define _DEBUG 37
#include <CGAL/Nef_3/debug.h>

CGAL_BEGIN_NAMESPACE

// ----------------------------------------------------------------------------
// SNC_ray_shooting
// ----------------------------------------------------------------------------

/*{\Manpage{SNC_ray_shooting}{SNC}{ray shoot functionality}{O}}*/

template <typename SNC_structure_>
class SNC_ray_shooter : public SNC_decorator<SNC_structure_>
{ 
  typedef SNC_structure_ SNC_structure;

protected:
  typedef SNC_ray_shooter<SNC_structure>          Self;
  typedef SNC_decorator<SNC_structure>            Base;

public:
  typedef typename SNC_structure_::Kernel          Kernel;
  typedef SNC_decorator<SNC_structure>             SNC_decorator;
  typedef SM_decorator<SNC_structure>          SM_decorator;
  typedef SM_point_locator<SM_decorator>       SM_point_locator;
  typedef SM_const_decorator<SNC_structure>    SM_const_decorator;
  typedef SNC_intersection<SNC_structure>          SNC_intersection;

  typedef typename SNC_structure::Vertex Vertex;
  typedef typename SNC_structure::Halfedge Halfedge;
  typedef typename SNC_structure::Halffacet Halffacet;
  typedef typename SNC_structure::Volume Volume;
  
  typedef typename SNC_structure::Vertex_iterator Vertex_iterator;
  typedef typename SNC_structure::Halfedge_iterator Halfedge_iterator;
  typedef typename SNC_structure::Halffacet_iterator Halffacet_iterator;
  typedef typename SNC_structure::Volume_iterator Volume_iterator;

  typedef typename SNC_structure::Vertex_handle Vertex_handle;
  typedef typename SNC_structure::Halfedge_handle Halfedge_handle;
  typedef typename SNC_structure::Halffacet_handle Halffacet_handle;
  typedef typename SNC_structure::Volume_handle Volume_handle;

  typedef typename SNC_structure::Vertex_const_handle Vertex_const_handle;
  typedef typename SNC_structure::Halfedge_const_handle Halfedge_const_handle;
  typedef typename SNC_structure::Halffacet_const_handle Halffacet_const_handle;
  typedef typename SNC_structure::Volume_const_handle Volume_const_handle;

  typedef typename SNC_structure::SVertex_iterator SVertex_iterator;
  typedef typename SNC_structure::SHalfedge_iterator SHalfedge_iterator;
  typedef typename SNC_structure::SFace_iterator SFace_iterator;
  typedef typename SNC_structure::SHalfloop_iterator SHalfloop_iterator;

  typedef typename SNC_structure::SVertex SVertex;
  typedef typename SNC_structure::SHalfedge SHalfedge;
  typedef typename SNC_structure::SFace SFace;
  typedef typename SNC_structure::SHalfloop SHalfloop;

  typedef typename SNC_structure::SVertex_handle SVertex_handle;
  typedef typename SNC_structure::SHalfedge_handle SHalfedge_handle;
  typedef typename SNC_structure::SFace_handle SFace_handle;
  typedef typename SNC_structure::SHalfloop_handle SHalfloop_handle;

  typedef typename SNC_structure::SVertex_const_handle SVertex_const_handle; 
  typedef typename SNC_structure::SHalfedge_const_handle SHalfedge_const_handle; 
  typedef typename SNC_structure::SHalfloop_const_handle SHalfloop_const_handle; 
  typedef typename SNC_structure::SFace_const_handle SFace_const_handle; 

  typedef typename SNC_structure::Object_handle Object_handle;
  typedef typename SNC_structure::SObject_handle SObject_handle;

  typedef typename SNC_structure::SHalfedge_around_facet_const_circulator SHalfedge_around_facet_const_circulator;
  typedef typename SNC_structure::SHalfedge_around_facet_circulator SHalfedge_around_facet_circulator;
  typedef typename SNC_structure::SFace_cycle_iterator SFace_cycle_iterator;
  typedef typename SNC_structure::SFace_cycle_const_iterator SFace_cycle_const_iterator;
  typedef typename SNC_structure::Halffacet_cycle_iterator Halffacet_cycle_iterator;
  typedef typename SNC_structure::Halffacet_cycle_const_iterator Halffacet_cycle_const_iterator;
  typedef typename SNC_structure::Shell_entry_iterator Shell_entry_iterator;
  typedef typename SNC_structure::Shell_entry_const_iterator Shell_entry_const_iterator;

  typedef typename SNC_structure::Point_3 Point_3;
  typedef typename SNC_structure::Vector_3 Vector_3;
  typedef typename SNC_structure::Segment_3 Segment_3;
  typedef typename SNC_structure::Ray_3 Ray_3;
  typedef typename SNC_structure::Line_3 Line_3;
  typedef typename SNC_structure::Plane_3 Plane_3;

  typedef typename SNC_structure::Sphere_point Sphere_point;
  typedef typename SNC_structure::Sphere_segment Sphere_segment;
  typedef typename SNC_structure::Sphere_circle Sphere_circle;
  typedef typename SNC_structure::Sphere_direction Sphere_direction;

  typedef typename SNC_structure::Mark Mark;
  typedef typename SNC_structure::Infi_box Infi_box;


  typedef typename SM_decorator::SHalfedge_around_svertex_const_circulator 
                                 SHalfedge_around_svertex_const_circulator;
  typedef typename SM_decorator::SHalfedge_around_svertex_circulator 
                                 SHalfedge_around_svertex_circulator;

  typedef void* GenPtr;

  SNC_ray_shooter() {}
  void initialize(SNC_structure* W) { Base::initialize(W); }

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
      TRACEN("facet below from from vertex...");
      f_below = get_visible_facet(v, ray);
      if(f_below != Halffacet_handle())
	return volume(f_below);
      SM_decorator SD(v);
      CGAL_assertion( SD.number_of_sfaces() == 1);
      return volume(SD.sfaces_begin());
    }
    else if( CGAL::assign(e, o)) {
      TRACEN("facet below from from edge...");
      f_below = get_visible_facet(e, ray);
      if(f_below != Halffacet_handle())
	return volume(f_below);
      SM_decorator SD(source(e));
      CGAL_assertion(SD.is_isolated(e));
      return volume(sface(e));
    }
    else if( CGAL::assign(f, o)) {
      TRACEN("facet below from from facet...");
      f_below = get_visible_facet(f, ray);
      CGAL_assertion( f_below != Halffacet_handle());
      return volume(f_below);
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

    TRACEN( "Shooting ray " << ray);
    Object_handle o;
    Vertex_handle v;
    CGAL_forall_vertices( v, *this->sncp()) {
      if ( ray.source() != point(v) && ray.has_on(point(v))) {
        if(hit && !Segment_3(ray.source(), end_of_seg).has_on(point(v)))
          continue;
        TRACEN("ray hit vertex case "<<point(v));
        end_of_seg = point(v);
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
          TRACEN("ray hit edge case " << segment(e) << " in " << q);
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
        TRACEN("ray hit facet "<<plane(f)<<" on "<<q);
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

    TRACEN( "Point locator for " << p);
    Vertex_handle v;
    CGAL_forall_vertices( v, *this->sncp()) {
      TRACEN("test vertex " << point(v));
      if ( p == point(v)) {
	TRACEN("on vertex.");
	return Object_handle(v);
      }
    }

    Halfedge_handle e;
    CGAL_forall_edges( e, *this->sncp()) {
      if ( is.does_contain_internally( segment(e), p) ) {
	TRACEN("on edge.");
	return Object_handle(e);
      }
    }
    Halffacet_handle f;
    CGAL_forall_halffacets( f, *this->sncp()) {
      if ( is.does_contain_internally( f, p) ) {
	TRACEN("on facet.");
	return Object_handle(f);
      }
    }

    Ray_3 r( p, Vector_3( -1, 0, 0));
    return Object_handle(determine_volume(r));
  }   

}; // SNC_ray_shooter

CGAL_END_NAMESPACE

#endif //CGAL_SNC_RAY_SHOOTER_H


