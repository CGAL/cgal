// ============================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : include/CGAL/Nef_3/SNC_ray_shooter.h
// package       : Nef_3
// chapter       : 3D-Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// maintainer    : Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// coordinator   : MPI Saarbruecken
//
// SNC_ray_shooter.h       Ray shooting and  point location on SNC structures
// ============================================================================
#ifndef CGAL_SNC_RAY_SHOOTER_H
#define CGAL_SNC_RAY_SHOOTER_H

#include <CGAL/basic.h>
#include <CGAL/functional.h> 
#include <CGAL/function_objects.h> 
#include <CGAL/Circulator_project.h>
#include <CGAL/Nef_3/bounded_side_3.h>
#include <CGAL/Nef_3/Pluecker_line_3.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_SM_overlayer.h>
#include <CGAL/Nef_3/SNC_SM_point_locator.h>
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
  typedef SNC_SM_decorator<SNC_structure>          SM_decorator;
  typedef SNC_SM_point_locator<SNC_structure>      SM_point_locator;
  typedef SNC_SM_const_decorator<SNC_structure>    SM_const_decorator;
  typedef SNC_intersection<SNC_structure>          SNC_intersection;

  #define USING(t) typedef typename SNC_structure::t t
  USING(Vertex);
  USING(Halfedge);
  USING(Halffacet);
  USING(Volume);
  
  USING(Vertex_iterator);
  USING(Halfedge_iterator);
  USING(Halffacet_iterator);
  USING(Volume_iterator);

  USING(Vertex_handle);
  USING(Halfedge_handle);
  USING(Halffacet_handle);
  USING(Volume_handle);

  USING(Vertex_const_handle);
  USING(Halfedge_const_handle);
  USING(Halffacet_const_handle);
  USING(Volume_const_handle);

  USING(SVertex_iterator);
  USING(SHalfedge_iterator);
  USING(SFace_iterator);
  USING(SHalfloop_iterator);

  USING(SVertex);
  USING(SHalfedge);
  USING(SFace);
  USING(SHalfloop);

  USING(SVertex_handle);
  USING(SHalfedge_handle);
  USING(SFace_handle);
  USING(SHalfloop_handle);

  USING(SVertex_const_handle); 
  USING(SHalfedge_const_handle); 
  USING(SHalfloop_const_handle); 
  USING(SFace_const_handle); 

  USING(Object_handle);
  USING(SObject_handle);

  USING(SHalfedge_around_facet_const_circulator);
  USING(SHalfedge_around_facet_circulator);
  USING(SFace_cycle_iterator);
  USING(SFace_cycle_const_iterator);
  USING(Halffacet_cycle_iterator);
  USING(Halffacet_cycle_const_iterator);
  USING(Shell_entry_iterator);
  USING(Shell_entry_const_iterator);

  USING(Point_3);
  USING(Vector_3);
  USING(Direction_3);
  USING(Segment_3);
  USING(Ray_3);
  USING(Line_3);
  USING(Plane_3);

  USING(Sphere_point);
  USING(Sphere_segment);
  USING(Sphere_circle);
  USING(Sphere_direction);

  USING(Mark);
  USING(Infi_box);
  #undef USING

  #define DECUSING(t) typedef typename SM_decorator::t t
  DECUSING(SHalfedge_around_svertex_const_circulator);
  DECUSING(SHalfedge_around_svertex_circulator);
  #undef DECUSING

  typedef void* GenPtr;

  SNC_ray_shooter(SNC_structure& W) : Base(W) {}
  /*{\Mcreate makes |\Mvar| a ray shooter on |W|.}*/

  Volume_handle determine_volume(Ray_3& ray) const {
    
    Object_handle o = shoot(ray);
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f, f_below;
    TRACEN("get_facet_below");
    if( assign(v, o)) {
      TRACEN("facet below from from vertex...");
      f_below = get_visible_facet(v, ray);
      if(f_below != Halffacet_handle())
	return volume(f_below);
      SM_decorator SD(v);
      CGAL_nef3_assertion( !is_empty_range( SD.sfaces_begin(), SD.sfaces_end()));
      CGAL_nef3_assertion( is_empty_range( ++(SD.sfaces_begin()), SD.sfaces_end()));
      return volume(SD.sfaces_begin());
    }
    else if( assign(e, o)) {
      TRACEN("facet below from from edge...");
      f_below = get_visible_facet(e, ray);
      if(f_below != Halffacet_handle())
	return volume(f_below);
      SM_decorator SD(source(e));
      CGAL_nef3_assertion(SD.is_isolated(e));
      return volume(sface(e));
    }
    else if( assign(f, o)) {
      TRACEN("facet below from from facet...");
      f_below = get_visible_facet(f, ray);
      CGAL_nef3_assertion( f_below != Halffacet_handle());
      return volume(f_below);
    }
    
    return Base(*this).volumes_begin();
  }

  Object_handle shoot( Ray_3& ray) const 
     /*{\Mop returns the nearest object hit by a ray |ray|. }*/ {
    bool hit = false;
    Point_3 end_of_seg;
    SNC_intersection is(*sncp());

    TRACEN( "Shooting ray " << ray);
    Object_handle o;
    Vertex_handle v;
    CGAL_nef3_forall_vertices( v, *sncp()) {
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
    CGAL_nef3_forall_edges( e, *sncp()) {
      Point_3 q;
      if( is.does_intersect_internally( ray, segment(e), q)) {
	if (!hit || has_smaller_distance_to_point(ray.source(),q, end_of_seg)) {
	  TRACEN("ray hit edge case " << segment(e) << " in " << q);
	  end_of_seg = q; 
	  hit = true;
	  o = Object_handle(e);
	}
      }
    }

    Halffacet_handle f;
    CGAL_nef3_forall_halffacets( f, *sncp()) {
      Point_3 q;
      if( is.does_intersect_internally( ray, f, q) ) {
	if(!hit || has_smaller_distance_to_point(ray.source(), q, end_of_seg)) {
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
    //{\Mop returns the lowest dimension object on an SNC structure
    //  which contais |p| in its interior. }
    {
      //            SETDTHREAD(37*19*47);

    SNC_intersection is(*sncp());

    TRACEN( "Point locator for " << p);
    Vertex_handle v;
    CGAL_nef3_forall_vertices( v, *sncp()) {
      TRACEN("test vertex " << point(v));
      if ( p == point(v)) {
	TRACEN("on vertex.");
	return Object_handle(v);
      }
    }

    Halfedge_handle e;
    CGAL_nef3_forall_edges( e, *sncp()) {
      if ( is.does_contain_internally( segment(e), p) ) {
	TRACEN("on edge.");
	return Object_handle(e);
      }
    }
    Halffacet_handle f;
    CGAL_nef3_forall_halffacets( f, *sncp()) {
      if ( is.does_contain_internally( f, p) ) {
	TRACEN("on facet.");
	return Object_handle(f);
      }
    }

    // let |s| be the segment that connects |p| to any fixed vertex |va| 
    Vertex_handle va = --(sncp()->vertices_end()); 
    TRACEN("target of shooting is " << point(va));
    Ray_3 s( p, point(va));
    return Object_handle(determine_volume(s));

    /*
    // SNC_decorator D(*this);
    //    Ray_3 s = Infi_box::get_ray(D ,p);
    // prune |s| by |o| if |o| intersects |s| in its relative interior 
    Object_handle o = shoot(s);
    // determine the volume that contains |s| from the last pruning object 

    Halffacet_handle hf;
    if( assign( v, o)) {
      hf = get_visible_facet( v, s);
      CGAL_nef3_assertion(hf != Halffacet_handle());
      //      if(hf==Halffacet_handle()) {
      //	SM_decorator SD(v);
      //	CGAL_nef3_assertion( !is_empty_range( SD.sfaces_begin(), SD.sfaces_end()));
      //	CGAL_nef3_assertion( is_empty_range( ++(SD.sfaces_begin()), SD.sfaces_end()));
      //	return volume(SD.sfaces_begin());
      //      }
    }
    else if( assign( e, o))
      hf = get_visible_facet( e, s);
    else if( assign( f, o)) {
      hf = get_visible_facet( f, s);
      CGAL_nef3_assertion(hf!=Halffacet_handle());
    }
    else CGAL_nef3_assertion_msg(0, "where is our point, eh?");
    
    if(hf!=Halffacet_handle())
      return Object_handle(volume(hf));

    if( assign( e, o))
      v = source(e);

    SM_decorator SD(v);
    CGAL_nef3_assertion( !is_empty_range( SD.sfaces_begin(), SD.sfaces_end()));
    CGAL_nef3_assertion( is_empty_range( ++(SD.sfaces_begin()), SD.sfaces_end()));
    return volume(SD.sfaces_begin());
    */
  }  

}; // SNC_ray_shoot

CGAL_END_NAMESPACE

#endif //CGAL_SNC_RAY_SHOOTER_H


