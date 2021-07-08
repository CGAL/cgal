// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Miguel Granados <granados@mpi-sb.mpg.de>

#ifndef CGAL_NEF_SNC_POINT_LOCATOR_H
#define CGAL_NEF_SNC_POINT_LOCATOR_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/basic.h>
#include <CGAL/Nef_3/SNC_intersection.h>
#include <CGAL/Nef_3/SNC_k3_tree_traits.h>
#include <CGAL/Nef_3/K3_tree.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Timer.h>
#include <string>


#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 509
#include <CGAL/Nef_2/debug.h>

#undef _CGAL_NEF_TRACEN
#define _CGAL_NEF_TRACEN(msg) CGAL_NEF_TRACEN( "SNC_point_locator: " << msg)

// TODO: find out the proper CGAL replacement for this macro and remove it
#define CGAL_for_each( i, C) for( i = C.begin(); i != C.end(); ++i)

// #define CGAL_NEF_TIMER(instruction) instruction
#define CGAL_NEF_TIMER(instruction)

// #define CGAL_NEF_CLOG(t) std::clog <<" "<<t<<std::endl; std::clog.flush()
#define CGAL_NEF_CLOG(t)

namespace CGAL {

template <typename SNC_decorator>
class SNC_point_locator
{
 public:
  class Intersection_call_back;
  typedef SNC_decorator Base;
  typedef SNC_point_locator<SNC_decorator> Self;
  typedef typename SNC_decorator::Decorator_traits Decorator_traits;
  typedef typename SNC_decorator::SNC_structure SNC_structure;
protected:
  std::string version_;
  // time for construction, point location, ray shooting and intersection test
  mutable Timer ct_t, pl_t, rs_t, it_t;

public:
  typedef typename SNC_structure::Object_handle Object_handle;

  typedef typename SNC_structure::Point_3 Point_3;
  typedef typename SNC_structure::Segment_3 Segment_3;
  typedef typename SNC_structure::Ray_3 Ray_3;
  typedef typename SNC_structure::Vector_3 Vector_3;
  typedef typename SNC_structure::Aff_transformation_3
                                  Aff_transformation_3;

  typedef typename Decorator_traits::Vertex_handle Vertex_handle;
  typedef typename Decorator_traits::Halfedge_handle Halfedge_handle;
  typedef typename Decorator_traits::Halffacet_handle Halffacet_handle;
  typedef typename Decorator_traits::Volume_handle Volume_handle;
  typedef typename Decorator_traits::Vertex_iterator Vertex_iterator;
  typedef typename Decorator_traits::Halfedge_iterator Halfedge_iterator;
  typedef typename Decorator_traits::Halffacet_iterator Halffacet_iterator;


  const std::string& version() const { return version_; }

  virtual Object_handle locate(const Point_3& p) const = 0;

  virtual Object_handle shoot(const Ray_3& s, int mask=255) const = 0;

  virtual void intersect_with_edges( Halfedge_handle edge,
                                     const Intersection_call_back& call_back)
    const = 0;

  virtual void intersect_with_facets( Halfedge_handle edge,
                                      const Intersection_call_back& call_back)
    const = 0;

  virtual void intersect_with_edges_and_facets( Halfedge_handle edge,
        const Intersection_call_back& call_back) const = 0;

  class Intersection_call_back
  {
  public:
    virtual void operator()( Halfedge_handle edge, Object_handle object,
                             const Point_3& intersection_point) const = 0;

    virtual ~Intersection_call_back() {}
  };

  virtual void initialize(SNC_structure* W) = 0;

  virtual Self* clone() const = 0;

  virtual void transform(const Aff_transformation_3& t) = 0;

  //virtual bool update( Unique_hash_map<Vertex_handle, bool>& V,
  //                   Unique_hash_map<Halfedge_handle, bool>& E,
  //                   Unique_hash_map<Halffacet_handle, bool>& F) = 0;

  virtual void add_facet(Halffacet_handle) {}

  virtual void add_edge(Halfedge_handle) {}

  virtual void add_vertex(Vertex_handle) {}

  virtual ~SNC_point_locator() noexcept(!CGAL_ASSERTIONS_ENABLED)
  {
    CGAL_NEF_CLOG("");
    CGAL_NEF_CLOG("construction_time:  "<<ct_t.time());
    CGAL_NEF_CLOG("pointlocation_time: "<<pl_t.time());
    CGAL_NEF_CLOG("rayshooting_time:   "<<rs_t.time());
    CGAL_NEF_CLOG("intersection_time:  "<<it_t.time());
    // warning: the total time showed here could be actually larger
    // that the real time used by this class, since point location
    // and intersection test use the ray shooter and so the same time
    // could be account to more than one timer
    CGAL_NEF_CLOG("pltotal_time:       "<<
      ct_t.time()+pl_t.time()+rs_t.time()+it_t.time());
  };
};

template <typename SNC_decorator>
class SNC_point_locator_by_spatial_subdivision :
  public SNC_point_locator<SNC_decorator>,
  public SNC_decorator
{
  typedef SNC_decorator Base;
  typedef CGAL::SNC_point_locator<SNC_decorator> SNC_point_locator;
  typedef CGAL::SNC_point_locator_by_spatial_subdivision<SNC_decorator> Self;
  typedef typename SNC_decorator::SNC_structure SNC_structure;
  typedef typename SNC_decorator::Decorator_traits Decorator_traits;
  typedef typename Decorator_traits::SM_decorator SM_decorator;
  typedef CGAL::SNC_intersection<SNC_structure> SNC_intersection;
  typedef typename SNC_decorator::SM_point_locator SM_point_locator;
  typedef typename SNC_decorator::Kernel Kernel;

public:
  typedef typename CGAL::SNC_k3_tree_traits<SNC_decorator> K3_tree_traits;

  typedef typename CGAL::K3_tree<K3_tree_traits> K3_tree;
  typedef K3_tree SNC_candidate_provider;

  typedef typename SNC_structure::Object_handle Object_handle;
  typedef typename SNC_structure::Point_3 Point_3;
  typedef typename SNC_structure::Plane_3 Plane_3;
  typedef typename SNC_structure::Segment_3 Segment_3;
  typedef typename SNC_structure::Ray_3 Ray_3;
  typedef typename SNC_structure::Vector_3 Vector_3;
  typedef typename SNC_structure::Triangle_3 Triangle_3;
  typedef typename SNC_structure::Aff_transformation_3
                                  Aff_transformation_3;

  typedef typename SNC_structure::Infi_box Infi_box;

  typedef typename Decorator_traits::Vertex_handle Vertex_handle;
  typedef typename Decorator_traits::Halfedge_handle Halfedge_handle;
  typedef typename Decorator_traits::Halffacet_handle Halffacet_handle;
  typedef typename Decorator_traits::SHalfedge_handle SHalfedge_handle;
  typedef typename Decorator_traits::SFace_handle SFace_handle;
  typedef typename Decorator_traits::Vertex_iterator Vertex_iterator;
  typedef typename Decorator_traits::Halfedge_iterator Halfedge_iterator;
  typedef typename Decorator_traits::Halffacet_iterator Halffacet_iterator;
  typedef typename Decorator_traits::Volume_handle Volume_handle;

  typedef typename Decorator_traits::Halffacet_cycle_iterator
                                     Halffacet_cycle_iterator;
  typedef typename Decorator_traits::SHalfedge_around_facet_circulator
                                     SHalfedge_around_facet_circulator;

  typedef typename SNC_candidate_provider::Object_list Object_list;
  typedef typename Object_list::iterator Object_list_iterator;
  typedef typename SNC_candidate_provider::Objects_along_ray Objects_along_ray;
  typedef typename Objects_along_ray::Iterator Objects_along_ray_iterator;

  //  typedef typename SNC_candidate_provider::Objects_around_box Objects_around_box;

  using Base::get_visible_facet;
public:
  SNC_point_locator_by_spatial_subdivision() :
    initialized(false), candidate_provider(0) {}


  virtual void initialize(SNC_structure* W) {
#ifdef CGAL_NEF_LIST_OF_TRIANGLES
    this->set_snc(*W);
    candidate_provider = new SNC_candidate_provider(W);
#else // CGAL_NEF_LIST_OF_TRIANGLES
    CGAL_NEF_TIMER(ct_t.start());
    this->version_ = std::string("Point Locator by Spatial Subdivision (tm)");
    CGAL_NEF_CLOG(version());
    CGAL_assertion( W != nullptr);
//    (Base) *this = SNC_decorator(*W);
    this->set_snc(*W);
    Object_list objects;
    Vertex_iterator v;
    Halfedge_iterator e;
    Halffacet_iterator f;
    CGAL_forall_vertices( v, *this->sncp())
      objects.push_back(make_object(Vertex_handle(v)));
    typename Object_list::difference_type v_end = objects.size();
    CGAL_forall_edges( e, *this->sncp())
      objects.push_back(make_object(Halfedge_handle(e)));
    CGAL_forall_facets( f, *this->sncp()) {
      objects.push_back(make_object(Halffacet_handle(f)));
    }
    if(initialized)
      delete candidate_provider;

    Object_list_iterator oli=objects.begin()+v_end;
    candidate_provider = new SNC_candidate_provider(objects,oli);

    // CGAL_NEF_TRACEN(*candidate_provider);
    CGAL_NEF_TIMER(ct_t.stop());
#endif // CGAL_NEF_LIST_OF_TRIANGLES
    initialized = true;
  }

  virtual Self* clone() const {
    return new Self;
  }

  virtual void transform(const Aff_transformation_3& t) {
    candidate_provider->transform(t);
  }

  virtual bool update( Unique_hash_map<Vertex_handle, bool>& V,
                       Unique_hash_map<Halfedge_handle, bool>& E,
                       Unique_hash_map<Halffacet_handle, bool>& F) {
    CGAL_NEF_TIMER(ct_t.start());
    CGAL_assertion( initialized);
    bool updated = candidate_provider->update( V, E, F);
    CGAL_NEF_TIMER(ct_t.stop());
    return updated;
  }

  virtual ~SNC_point_locator_by_spatial_subdivision() noexcept(!CGAL_ASSERTIONS_ENABLED)
  {
    CGAL_destructor_warning(initialized ||
                 candidate_provider == 0); // required?
    if(initialized)
      delete candidate_provider;
  }

  /*
  template <typename Box>
  bool point_in_box(const Box& box, const Point_3& p) {
    if(p.x() < box.min_coord(0)) return true;
    if(p.x() > box.max_coord(0)) return true;
    if(p.y() < box.min_coord(1)) return true;
    if(p.y() > box.max_coord(1)) return true;
    if(p.z() < box.min_coord(2)) return true;
    if(p.z() > box.max_coord(2)) return true;
    return false;
  }

  template <typename Box>
  bool in_range(const Box& box) const {
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;

    Objects_around_box oab;
    typename Objects_around_box::Iterator
      objects_iterator(oab.begin());
    for( ; objects_iterator != oab.end(); ++objects_iterator) {
      Object_list candidates = *objects_iterator;
      Object_list_iterator o;
      CGAL_for_each( o, candidates) {
        if(CGAL::assign(v, *o)) {
          if(point_in_box(box, v->point()))
            return true;
        } else if(CGAL::assign(e, *o)) {
          Point_3 ip;
          Point_3 pmin(box.min_coord(0), box.min_coord(1), box.min_coord(2));
          Point_3 pmax(box.max_coord(0), box.max_coord(1), box.max_coord(2));
          Segment_3 seg(e->source()->point(),
                        e->twin()->source()->point());

          Plane_3 pl(pmin, Vector_3(1,0,0));
          Object o = intersection(pl, seg);
          if(CGAL::assign(ip,o) && point_in_box(box, ip))
            return true;
          pl = Plane_3(pmin, Vector_3(0,1,0));
          o = intersection(pl, seg);
          if(CGAL::assign(ip,o) && point_in_box(box, ip))
            return true;
          pl = Plane_3(pmin, Vector_3(0,0,1));
          o = intersection(pl, seg);
          if(CGAL::assign(ip,o) && point_in_box(box, ip))
            return true;
          pl = Plane_3(pmax, Vector_3(1,0,0));
          o = intersection(pl, seg);
          if(CGAL::assign(ip,o) && point_in_box(box, ip))
            return true;
          pl = Plane_3(pmax, Vector_3(0,1,0));
          o = intersection(pl, seg);
          if(CGAL::assign(ip,o) && point_in_box(box, ip))
            return true;
          pl = Plane_3(pmax, Vector_3(0,0,1));
          o = intersection(pl, seg);
          if(CGAL::assign(ip,o) && point_in_box(box, ip))
            return true;
        } else if(CGAL::assign(f, *o)) {
          Segment_3 seg(Point_3(box.min_coord(0), box.min_coord(1), box.min_coord(2)),
                        Point_3(box.min_coord(0), box.min_coord(1), box.max_coord(2)));
          if(is.does_intersect_internally(seg, f))
            return true;
          seg = Segment_3(Point_3(box.max_coord(0), box.min_coord(1), box.min_coord(2)),
                          Point_3(box.max_coord(0), box.min_coord(1), box.max_coord(2)));
          if(is.does_intersect_internally(seg, f))
            return true;
          seg = Segment_3(Point_3(box.min_coord(0), box.max_coord(1), box.min_coord(2)),
                          Point_3(box.min_coord(0), box.max_coord(1), box.max_coord(2)));
          if(is.does_intersect_internally(seg, f))
            return true;
          seg = Segment_3(Point_3(box.max_coord(0), box.max_coord(1), box.min_coord(2)),
                          Point_3(box.max_coord(0), box.max_coord(1), box.max_coord(2)));
          if(is.does_intersect_internally(seg, f))
            return true;

          seg = Segment_3(Point_3(box.min_coord(0), box.min_coord(1), box.min_coord(2)),
                          Point_3(box.min_coord(0), box.max_coord(1), box.min_coord(2)));
          if(is.does_intersect_internally(seg, f))
            return true;
          seg = Segment_3(Point_3(box.max_coord(0), box.min_coord(1), box.min_coord(2)),
                          Point_3(box.max_coord(0), box.max_coord(1), box.min_coord(2)));
          if(is.does_intersect_internally(seg, f))
            return true;
          seg = Segment_3(Point_3(box.min_coord(0), box.min_coord(1), box.max_coord(2)),
                          Point_3(box.min_coord(0), box.max_coord(1), box.max_coord(2)));
          if(is.does_intersect_internally(seg, f))
            return true;
          seg = Segment_3(Point_3(box.max_coord(0), box.min_coord(1), box.max_coord(2)),
                          Point_3(box.max_coord(0), box.max_coord(1), box.max_coord(2)));
          if(is.does_intersect_internally(seg, f))
            return true;

          seg = Segment_3(Point_3(box.min_coord(0), box.min_coord(1), box.min_coord(2)),
                          Point_3(box.max_coord(0), box.min_coord(1), box.min_coord(2)));
          if(is.does_intersect_internally(seg, f))
            return true;
          seg = Segment_3(Point_3(box.min_coord(0), box.max_coord(1), box.min_coord(2)),
                          Point_3(box.max_coord(0), box.max_coord(1), box.min_coord(2)));
          if(is.does_intersect_internally(seg, f))
            return true;
          seg = Segment_3(Point_3(box.min_coord(0), box.min_coord(1), box.max_coord(2)),
                          Point_3(box.max_coord(0), box.min_coord(1), box.max_coord(2)));
          if(is.does_intersect_internally(seg, f))
            return true;
          seg = Segment_3(Point_3(box.min_coord(0), box.max_coord(1), box.max_coord(2)),
                          Point_3(box.max_coord(0), box.max_coord(1), box.max_coord(2)));
          if(is.does_intersect_internally(seg, f))
            return true;

        } else
          CGAL_assertion_msg(false, "wrong handle");
      }
    }
    return false;
  }
  */

  virtual Object_handle shoot(const Ray_3& ray, int mask=255) const {
    CGAL_NEF_TIMER(rs_t.start());
    CGAL_assertion( initialized);
    _CGAL_NEF_TRACEN( "shooting: "<<ray);
    Object_handle result;
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    bool hit = false;
    Point_3 eor = CGAL::ORIGIN; // 'end of ray', the latest ray's hit point
    Objects_along_ray objects = candidate_provider->objects_along_ray(ray);
    Objects_along_ray_iterator objects_iterator = objects.begin();
    while( !hit && objects_iterator != objects.end()) {
      Object_list candidates = *objects_iterator;
      Object_list_iterator o;
      CGAL_for_each( o, candidates) {
        if( CGAL::assign( v, *o) && ((mask&1) != 0)) {
          _CGAL_NEF_TRACEN("trying vertex on "<<v->point());
          if( ray.source() != v->point() && ray.has_on(v->point())) {
            _CGAL_NEF_TRACEN("the ray intersects the vertex");
            _CGAL_NEF_TRACEN("prev. intersection? "<<hit);
            CGAL_assertion_code
              (if( hit)_CGAL_NEF_TRACEN("prev. intersection on "<<eor));
            if( hit && !Segment_3( ray.source(), eor).has_on(v->point()))
              continue;
            eor = v->point();
            result = make_object(v);
            hit = true;
            _CGAL_NEF_TRACEN("the vertex becomes the new hit object");
          }
        }
        else if( CGAL::assign( e, *o) && ((mask&2) != 0)) {
          Point_3 q;
          _CGAL_NEF_TRACEN("trying edge on "<< Segment_3(e->source()->point(),e->twin()->source()->point()));
          if( is.does_intersect_internally( ray, Segment_3(e->source()->point(),
                                                           e->twin()->source()->point()), q)) {
            _CGAL_NEF_TRACEN("ray intersects edge on "<<q);
            _CGAL_NEF_TRACEN("prev. intersection? "<<hit);
            CGAL_assertion_code
              (if( hit) _CGAL_NEF_TRACEN("prev. intersection on "<<eor));
            if( hit && !has_smaller_distance_to_point( ray.source(), q, eor))
              continue;
            _CGAL_NEF_TRACEN("is the intersection point on the current cell? "<<
                    candidate_provider->is_point_on_cell( q, objects_iterator));
            if( !candidate_provider->is_point_on_cell( q, objects_iterator))
                continue;
            eor = q;
            result = make_object(e);
            hit = true;
            _CGAL_NEF_TRACEN("the edge becomes the new hit object");
          }
        }
        else if( CGAL::assign( f, *o) && ((mask&4) != 0)) {
          Point_3 q;
          _CGAL_NEF_TRACEN("trying facet with on plane "<<f->plane()<<
                  " with point on "<<f->plane().point());
          if( is.does_intersect_internally( ray, f, q, true) ) {
            _CGAL_NEF_TRACEN("ray intersects facet on "<<q);
            _CGAL_NEF_TRACEN("prev. intersection? "<<hit);
            if( hit) { _CGAL_NEF_TRACEN("prev. intersection on "<<eor); }
            if( hit && !has_smaller_distance_to_point( ray.source(), q, eor))
              continue;
            _CGAL_NEF_TRACEN("is the intersection point on the current cell? "<<
                    candidate_provider->is_point_on_cell( q, objects_iterator));
            if( !candidate_provider->is_point_on_cell( q, objects_iterator))
              continue;
            eor = q;
            result = make_object(f);
            hit = true;
            _CGAL_NEF_TRACEN("the facet becomes the new hit object");
          }
        }
        else if((mask&15) == 15)
          CGAL_error_msg( "wrong handle");
      }
      if(!hit)
        ++objects_iterator;
    }
    CGAL_NEF_TIMER(rs_t.stop());
    return result;
  }

  virtual Object_handle locate( const Point_3& p) const {
    if(Infi_box::extended_kernel()) {
    CGAL_NEF_TIMER(pl_t.start());
    CGAL_assertion( initialized);
    _CGAL_NEF_TRACEN( "locate "<<p);
    Object_handle result;
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    Object_list candidates = candidate_provider->objects_around_point(p);
    Object_list_iterator o = candidates.begin();
    bool found = false;
    while( !found && o != candidates.end()) {
      if( CGAL::assign( v, *o)) {
        if ( p == v->point()) {
          _CGAL_NEF_TRACEN("found on vertex "<<v->point());
          result = make_object(v);
          found = true;
        }
      }
      else if( CGAL::assign( e, *o)) {
        if ( is.does_contain_internally(Segment_3(e->source()->point(),e->twin()->source()->point()), p) ) {
          _CGAL_NEF_TRACEN("found on edge "<<Segment_3(e->source()->point(),e->twin()->source()->point()));
          result = make_object(e);
          found = true;
        }
      }
      else if( CGAL::assign( f, *o)) {
        if (is.does_contain_internally( f, p) ) {
          _CGAL_NEF_TRACEN("found on facet...");
          result = make_object(f);
          found = true;
        }
      }
      o++;
    }
    if( !found) {
      _CGAL_NEF_TRACEN("point not found in 2-skeleton");
      _CGAL_NEF_TRACEN("shooting ray to determine the volume");
      Ray_3 r( p, Vector_3( -1, 0, 0));
      result = make_object(determine_volume(r));
    }    CGAL_NEF_TIMER(pl_t.start());
    CGAL_NEF_TIMER(pl_t.stop());
    return result;


  } else {   // standard kernel


    CGAL_assertion( initialized);
    _CGAL_NEF_TRACEN( "locate "<<p);
    typename SNC_structure::FT min_distance;
    typename SNC_structure::FT tmp_distance;
    Object_handle result;
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    Object_list candidates = candidate_provider->objects_around_point(p);
    Object_list_iterator o = candidates.begin();

    if(candidates.empty())
      return make_object(Base(*this).volumes_begin());

    CGAL::assign(v,*o);
    CGAL_assertion(CGAL::assign(v,*o));
    if(p==v->point())
      return make_object(v);

    min_distance = CGAL::squared_distance(v->point(),p);
    result = make_object(v);
    ++o;
    while(o!=candidates.end() && CGAL::assign(v,*o)) {
      if ( p == v->point()) {
        _CGAL_NEF_TRACEN("found on vertex "<<v->point());
        return make_object(v);
      }
      tmp_distance = CGAL::squared_distance(v->point(),p);
      if(tmp_distance < min_distance) {
        result = make_object(v);
        min_distance = tmp_distance;
      }
      ++o;
    }

    CGAL::assign(v, result);
    Segment_3 s(p,v->point());
    // bool first = true;
    Point_3 ip;

    /*
    // TODO: das geht effizienter
    Object_list_iterator of(o);
    while(of != candidates.end() && assign(e, *of)) ++of;

    typename SNC_structure::SHalfedge_iterator sei;
    for(sei=v->shalfedges_begin(); sei!=v->shalfedges_end(); ++sei){
      if(sei->is_twin()) continue;
      Halffacet_handle fout = sei->facet();
      if(fout->is_twin()) fout = fout->twin();
      Object_list_iterator ofc(of);
      for(;ofc!=candidates.end();++ofc) {
        if(CGAL::assign(f,*ofc)) {
          if(f == fout->twin())
            std::cerr << "shit" << std::endl;
          if(f == fout) {
            Object_list_iterator oe(ofc);
            --ofc;
            candidates.erase(oe);
          }
        }
      }
    }
    */
    for(;o!=candidates.end();++o) {
      if( CGAL::assign( e, *o)) {
        //        if(first &&
        //           (e->source() == v  || e->twin()->source() == v)) continue;
        Segment_3 ss(e->source()->point(),e->twin()->source()->point());
        CGAL_NEF_TRACEN("test edge " << e->source()->point() << "->" << e->twin()->source()->point());
        if (is.does_contain_internally(ss, p)) {
        _CGAL_NEF_TRACEN("found on edge "<< ss);
          return make_object(e);
        }
        if((e->source() != v)  && (e->twin()->source() != v) && is.does_intersect_internally(s, ss, ip)) {
          // first = false;
          s = Segment_3(p, normalized(ip));
          result = make_object(e);
        }

      } else
      if( CGAL::assign( f, *o)) {
        CGAL_NEF_TRACEN("test facet " << f->plane());
        if (is.does_contain_internally(f,p) ) {
          _CGAL_NEF_TRACEN("found on facet...");
          return make_object(f);
        }

        // We next check if v is a vertex on the face to avoid a geometric test
        bool v_vertex_of_f = false;
        Halffacet_cycle_iterator fci;
        for(fci=f->facet_cycles_begin(); (! v_vertex_of_f) && (fci!=f->facet_cycles_end()); ++fci) {
          if(fci.is_shalfedge()) {
            SHalfedge_around_facet_circulator sfc(fci), send(sfc);
            CGAL_For_all(sfc,send) {
              if(sfc->source()->center_vertex() ==  v){
                v_vertex_of_f = true;
                break;
              }
            }
          }
        }


        if( (! v_vertex_of_f) &&  is.does_intersect_internally(s,f,ip) ) {
          s = Segment_3(p, normalized(ip));
          result = make_object(f);
        }
      }
      else CGAL_error_msg( "wrong handle type");
    }

    //CGAL_warning("altered code in SNC_point_locator");
    /*
      Halffacet_iterator fc;
      CGAL_forall_facets(fc, *this->sncp()) {
        CGAL_assertion(!is.does_intersect_internally(s,f,ip));
      }

      Halfedge_iterator ec;
      CGAL_forall_edges(ec, *this->sncp()) {
        Segment_3 ss(ec->source()->point(), ec->twin()->source()->point());
        CGAL_assertion(!is.does_intersect_internally(s,ss,ip));
      }

      Vertex_iterator vc;
      CGAL_forall_vertices(vc, *this->sncp()) {
        std::cerr << "test vertex " << vc->point() << std::endl;
        CGAL_assertion(vc->point() == s.target() || !s.has_on(vc->point()));
      }
    */

    if( CGAL::assign( v, result)) {
      _CGAL_NEF_TRACEN("vertex hit, obtaining volume..." << v->point());

      //CGAL_warning("altered code in SNC_point_locator");
      SM_point_locator L(&*v);
      //      Object_handle so = L.locate(s.source()-s.target(), true);
      Object_handle so = L.locate(s.source()-s.target());
      SFace_handle sf;
      if(CGAL::assign(sf,so))
        return make_object(sf->volume());
      CGAL_error_msg( "wrong handle type");
      return Object_handle();
/*
      SHalfedge_handle se;
      CGAL_assertion(CGAL::assign(se,so));
      CGAL_NEF_TRACEN("intersect segment " << s << " with edges");
      for(;ox!=candidates.end();++ox) {
        if(!CGAL::assign(e,*ox)) continue;
        CGAL_NEF_TRACEN("test edge " << e->source()->point() << "->" << e->twin()->source()->point());
        if(is.does_intersect_internally(s,Segment_3(e->source()->point(),e->twin()->source()->point()),ip)) {
          s = Segment_3(p, normalized(ip));
          result = make_object(e);
        }
      }
      CGAL_assertion(CGAL::assign(e,result));
      CGAL::assign(e,result);
      f = get_visible_facet(e, Ray_3(p, s.target()));
      if( f != Halffacet_handle())
        return f->incident_volume();
      SM_decorator SD(&*v); // now, the vertex has no incident facets
      CGAL_assertion( SD.number_of_sfaces() == 1);
      return SD.sfaces_begin()->volume();
*/
    } else if( CGAL::assign( f, result)) {
      _CGAL_NEF_TRACEN("facet hit, obtaining volume...");
      if(f->plane().oriented_side(p) == ON_NEGATIVE_SIDE)
        f = f->twin();
      return make_object(f->incident_volume());
    } else if( CGAL::assign(e, result)) {
      SM_decorator SD(&*e->source());
      if( SD.is_isolated(e))
        return make_object(e->incident_sface()->volume());
      return make_object(get_visible_facet(e,Ray_3(s.source(),s.to_vector()))->incident_volume());
    }
    CGAL_error_msg( "wrong handle type");
    return Object_handle();
  }
  }

  virtual void intersect_with_edges_and_facets( Halfedge_handle e0,
        const typename SNC_point_locator::Intersection_call_back& call_back) const {

    CGAL_NEF_TIMER(it_t.start());
    CGAL_assertion( initialized);
    _CGAL_NEF_TRACEN( "intersecting edge: "<<&*e0<<' '<<Segment_3(e0->source()->point(),
                                                         e0->twin()->source()->point()));


    Segment_3 s(Segment_3(e0->source()->point(),e0->twin()->source()->point()));
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    Object_list_iterator o;
    Object_list objects = candidate_provider->objects_around_segment(s);
    CGAL_for_each( o, objects) {
      if( CGAL::assign( v, *o)) {
        /* do nothing */
      }
      else if( CGAL::assign( e, *o)) {

#ifdef CGAL_NEF3_DUMP_STATISTICS
      ++number_of_intersection_candidates;
#endif

        Point_3 q;
        if( is.does_intersect_internally( s, Segment_3(e->source()->point(),
                                                       e->twin()->source()->point()), q)) {
          q = normalized(q);
          call_back( e0, make_object(Halfedge_handle(e)), q);
          _CGAL_NEF_TRACEN("edge intersects edge "<<' '<<&*e<< Segment_3(e->source()->point(),
                                                                e->twin()->source()->point())<<" on "<<q);
        }
      }
      else if( CGAL::assign( f, *o)) {
#ifdef CGAL_NEF3_DUMP_STATISTICS
      ++number_of_intersection_candidates;
#endif

        Point_3 q;
        if( is.does_intersect_internally( s, f, q) ) {
          q = normalized(q);
          call_back( e0, make_object(Halffacet_handle(f)), q);
          _CGAL_NEF_TRACEN("edge intersects facet on plane "<<f->plane()<<" on "<<q);
        }
      }
      else
        CGAL_error_msg( "wrong handle");
    }
    CGAL_NEF_TIMER(it_t.stop());
  }

  virtual void intersect_with_edges( Halfedge_handle e0,
    const typename SNC_point_locator::Intersection_call_back& call_back) const {
    CGAL_NEF_TIMER(it_t.start());
    CGAL_assertion( initialized);
    _CGAL_NEF_TRACEN( "intersecting edge: "<<&*e0<<' '<<Segment_3(e0->source()->point(),
                                                         e0->twin()->source()->point()));
    Segment_3 s(Segment_3(e0->source()->point(),e0->twin()->source()->point()));
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    Object_list_iterator o;
    Object_list objects = candidate_provider->objects_around_segment(s);
    CGAL_for_each( o, objects) {
      if( CGAL::assign( v, *o)) {
        /* do nothing */
      }
      else if( CGAL::assign( e, *o)) {

#ifdef CGAL_NEF3_DUMP_STATISTICS
      ++number_of_intersection_candidates;
#endif

        Point_3 q;
        if( is.does_intersect_internally( s, Segment_3(e->source()->point(),
                                                       e->twin()->source()->point()), q)) {
          q = normalized(q);
          call_back( e0, make_object(Halfedge_handle(e)), q);
          _CGAL_NEF_TRACEN("edge intersects edge "<<' '<<&*e<< Segment_3(e->source()->point(),
                                                                e->twin()->source()->point())<<" on "<<q);
        }
      }
      else if( CGAL::assign( f, *o)) {
        /* do nothing */
      }
      else
        CGAL_error_msg( "wrong handle");
    }
    CGAL_NEF_TIMER(it_t.stop());
  }

  virtual void intersect_with_facets( Halfedge_handle e0,
    const typename SNC_point_locator::Intersection_call_back& call_back) const {
    CGAL_NEF_TIMER(it_t.start());
    CGAL_assertion( initialized);
    _CGAL_NEF_TRACEN( "intersecting edge: "<< Segment_3(e0->source()->point(),
                                               e0->twin()->source()->point()));
    Segment_3 s(Segment_3(e0->source()->point(),e0->twin()->source()->point()));
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    Object_list_iterator o;
    Object_list objects = candidate_provider->objects_around_segment(s);
    CGAL_for_each( o, objects) {
      if( CGAL::assign( v, *o)) {
        /* do nothing */
      }
      else if( CGAL::assign( e, *o)) {
        /* do nothing */
      }
      else if( CGAL::assign( f, *o)) {

#ifdef CGAL_NEF3_DUMP_STATISTICS
      ++number_of_intersection_candidates;
#endif

        Point_3 q;
        if( is.does_intersect_internally( s, f, q) ) {
          q = normalized(q);
          call_back( e0, make_object(Halffacet_handle(f)), q);
          _CGAL_NEF_TRACEN("edge intersects facet on plane "<<f->plane()<<" on "<<q);
        }
      }
      else
        CGAL_error_msg( "wrong handle");
    }
    CGAL_NEF_TIMER(it_t.stop());
  }

private:
  Volume_handle determine_volume( const Ray_3& ray) const {
    Halffacet_handle f_below;
    Object_handle o = shoot(ray);
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    if( CGAL::assign( v, o)) {
      _CGAL_NEF_TRACEN("vertex hit, obtaining volume...");
      f_below = get_visible_facet( v, ray);
      if( f_below != Halffacet_handle())
        return f_below->incident_volume();
      SM_decorator SD(&*v); // now, the vertex has no incident facets
      CGAL_assertion( SD.number_of_sfaces() == 1);
      return SD.sfaces_begin()->volume();
    }
    else if( CGAL::assign( e, o)) {
      _CGAL_NEF_TRACEN("edge hit, obtaining volume...");
      f_below = get_visible_facet( e, ray);
      if( f_below != Halffacet_handle())
        return f_below->incident_volume();
      CGAL_assertion_code(SM_decorator SD(&*e->source())); // now, the edge has no incident facets
      CGAL_assertion(SD.is_isolated(e));
      return e->incident_sface()->volume();
    }
    else if( CGAL::assign( f, o)) {
      _CGAL_NEF_TRACEN("facet hit, obtaining volume...");
      f_below = get_visible_facet(f, ray);
      CGAL_assertion( f_below != Halffacet_handle());
      return f_below->incident_volume();
    }
    return Base(*this).volumes_begin(); // TODO: Comment this hack!
  }

public:
  void add_facet(Halffacet_handle f) {
    candidate_provider->add_facet(f);
  }

  void add_edge(Halfedge_handle e) {
    candidate_provider->add_edge(e);
  }

  void add_vertex(Vertex_handle v) {
    candidate_provider->add_vertex(v);
  }

private:
  bool initialized;
  SNC_candidate_provider* candidate_provider;
  SNC_intersection is;
};

#ifdef CGAL_NEF3_POINT_LOCATOR_NAIVE
template <typename SNC_decorator>
class SNC_point_locator_naive :
  public SNC_ray_shooter<SNC_decorator>,
  public SNC_point_locator<SNC_decorator>
{
  typedef typename SNC_decorator::SNC_structure SNC_structure;
  typedef SNC_ray_shooter<SNC_decorator> Base;
  typedef SNC_point_locator_naive<SNC_decorator> Self;
  typedef SNC_point_locator<SNC_decorator> SNC_point_locator;
  typedef CGAL::SNC_intersection<SNC_structure> SNC_intersection;
  typedef typename SNC_decorator::Decorator_traits Decorator_traits;
  typedef typename Decorator_traits::SM_decorator SM_decorator;

public:
  typedef typename SNC_decorator::Object_handle Object_handle;
  typedef typename SNC_decorator::Point_3 Point_3;
  typedef typename SNC_decorator::Segment_3 Segment_3;
  typedef typename SNC_decorator::Ray_3 Ray_3;
  typedef typename SNC_structure::Aff_transformation_3
                                  Aff_transformation_3;

  typedef typename Decorator_traits::Vertex_handle Vertex_handle;
  typedef typename Decorator_traits::Halfedge_handle Halfedge_handle;
  typedef typename Decorator_traits::Halffacet_handle Halffacet_handle;
  typedef typename Decorator_traits::Volume_handle Volume_handle;
  typedef typename Decorator_traits::Vertex_iterator Vertex_iterator;
  typedef typename Decorator_traits::Halfedge_iterator Halfedge_iterator;
  typedef typename Decorator_traits::Halffacet_iterator Halffacet_iterator;

public:
  SNC_point_locator_naive() : initialized(false) {}
  virtual void initialize(SNC_structure* W) {
    CGAL_NEF_TIMER(ct_t.start());
    this->version_ = std::string("Naive Point Locator (tm)");
    CGAL_NEF_CLOG(version());
    CGAL_assertion( W != nullptr);
    Base::initialize(W);
    initialized = true;
    CGAL_NEF_TIMER(ct_t.stop());
  }

  virtual Self* clone() const {
    return new Self;
  }

  virtual bool update( Unique_hash_map<Vertex_handle, bool>& V,
                       Unique_hash_map<Halfedge_handle, bool>& E,
                       Unique_hash_map<Halffacet_handle, bool>& F) {
    CGAL_NEF_TIMER(ct_t.start());
    CGAL_assertion( initialized);
    CGAL_NEF_TIMER(ct_t.stop());
    return false;
  }

  virtual ~SNC_point_locator_naive() {}

  virtual Object_handle locate(const Point_3& p) const {
    CGAL_NEF_TIMER(pl_t.start());
    CGAL_assertion( initialized);
    CGAL_NEF_TIMER(pl_t.stop());
    return Base::locate(p);
  }

  virtual Object_handle shoot(const Ray_3& r, int mask=0) const {
    CGAL_NEF_TIMER(rs_t.start());
    CGAL_assertion( initialized);
    CGAL_NEF_TIMER(rs_t.stop());
    return Base::shoot(r);
  }

  virtual void transform(const Aff_transformation_3& aff) {}

  virtual void intersect_with_edges_and_facets(Halfedge_handle e0,
    const typename SNC_point_locator::Intersection_call_back& call_back) const {
        intersect_with_edges(e0,call_back);
        intersect_with_facets(e0,call_back);
  }

  virtual void intersect_with_edges( Halfedge_handle e0,
    const typename SNC_point_locator::Intersection_call_back& call_back) const {
    CGAL_NEF_TIMER(it_t.start());
    CGAL_assertion( initialized);
    CGAL_NEF_TRACEN( "intersecting edge: "<< Segment_3(e0->source()->point(),
                                              e0->twin()->source()->point()));
    SNC_intersection is(*this->sncp());
    Segment_3 s(Segment_3(e0->source()->point(),e0->twin()->source()->point()));
    Halfedge_iterator e;
    CGAL_forall_edges( e, *this->sncp()) {

#ifdef CGAL_NEF3_DUMP_STATISTICS
      ++number_of_intersection_candidates;
#endif

      Point_3 q;
      if( is.does_intersect_internally( s, Segment_3(e->source()->point(),
                                                     e->twin()->source()->point()), q)) {
        q = normalized(q);
        CGAL_NEF_TRACEN("edge intersects edge "<< Segment_3(e->source()->point(),
                                                   e->twin()->source()->point()) <<" on "<<q);
        call_back( e0, make_object(e), q);
      }
    }
    CGAL_NEF_TIMER(it_t.stop());
  }

  virtual void intersect_with_facets( Halfedge_handle e0,
    const typename SNC_point_locator::Intersection_call_back& call_back) const {
    CGAL_NEF_TIMER(it_t.start());
    CGAL_assertion( initialized);
    CGAL_NEF_TRACEN( "intersecting edge: "<< Segment_3(e0->source()->point(),
                                              e0->twin()->source()->point()));
    SNC_intersection is(*this->sncp());
    Segment_3 s(Segment_3(e0->source()->point(),
                          e0->twin()->source()->point()));
    Halffacet_iterator f;
    CGAL_forall_facets( f, *this->sncp()) {

#ifdef CGAL_NEF3_DUMP_STATISTICS
      ++number_of_intersection_candidates;
#endif

      Point_3 q;
      if( is.does_intersect_internally( s, f, q) ) {
        q = normalized(q);
        CGAL_NEF_TRACEN("edge intersects facet on plane "<<f->plane()<<" on "<<q);
        call_back( e0, make_object(f), q);
      }
    }
    CGAL_NEF_TIMER(it_t.stop());
  }

private:
  bool initialized;
};
#endif

} //namespace CGAL
#endif // CGAL_NEF_SNC_POINT_LOCATOR_H
