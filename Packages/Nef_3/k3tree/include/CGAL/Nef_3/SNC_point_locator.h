#line 7 "point_locator.nw"
#ifndef SNC_POINT_LOCATOR_H
#define SNC_POINT_LOCATOR_H

#include <CGAL/basic.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_intersection.h>
#include <CGAL/Nef_3/SNC_ray_shooter.h>
#include <CGAL/Nef_3/SNC_k3_tree_traits.h>
#include <CGAL/Nef_3/K3_tree.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Timer.h>

// #include <CGAL/Polygon_triangulation_traits_2.h>
// #include <CGAL/Nef_3/triangulate_nef3_facet.h>

#undef _DEBUG
#define _DEBUG 509
#include <CGAL/Nef_3/debug.h>

#undef _TRACEN
#define _TRACEN(msg) TRACEN( "SNC_point_locator: " << msg);

// TODO: find out the proper CGAL replacement for this macro and remove it
#define CGAL_for_each( i, C) for( i = C.begin(); i != C.end(); ++i)

// #define TIMER(instruction) instruction
#define TIMER(instruction)

// #define CLOG(t) std::clog <<" "<<t<<std::endl; std::clog.flush()
#define CLOG(t)

CGAL_BEGIN_NAMESPACE

template <typename SNC_structure>
class SNC_point_locator
{
 public:
  class Intersection_call_back;
  typedef SNC_point_locator<SNC_structure> Self;

protected:
  char version_[64];
  // time for construction, point location, ray shooting and intersection test
  mutable Timer ct_t, pl_t, rs_t, it_t; 

public: 
  #define USING(t) typedef typename SNC_structure::t t
  USING(Object_handle);
  USING(Vertex_handle);
  USING(Halfedge_handle);
  USING(Halffacet_handle);
  USING(Halffacet_triangle_handle);
  USING(Volume_handle);
  USING(Vertex_iterator);
  USING(Halfedge_iterator);
  USING(Halffacet_iterator);
  USING(Point_3);
  USING(Segment_3);
  USING(Ray_3);
  USING(Vector_3);
  USING(Triangle_3);
  USING(Aff_transformation_3);
  #undef USING

  const char* version() { return version_; }

  virtual Object_handle locate(const Point_3& p) const = 0;

  virtual Object_handle shoot(const Ray_3& s) const = 0;

  virtual void intersect_with_edges( Halfedge_handle edge,
                                     const Intersection_call_back& call_back) 
    const = 0;

  virtual void intersect_with_facets( Halfedge_handle edge,
                                      const Intersection_call_back& call_back)
    const = 0;

  class Intersection_call_back 
  {
  public:
    virtual void operator()( Halfedge_handle edge, Object_handle object, 
                             const Point_3& intersection_point) const = 0;
  };

  virtual void initialize(SNC_structure* W) = 0;

  virtual Self* clone() const = 0;

  virtual void transform(const Aff_transformation_3& t) = 0;

  //virtual bool update( Unique_hash_map<Vertex_handle, bool>& V, 
  //                   Unique_hash_map<Halfedge_handle, bool>& E, 
  //                   Unique_hash_map<Halffacet_handle, bool>& F) = 0;

  virtual ~SNC_point_locator() {
    CLOG("");
    CLOG("construction_time:  "<<ct_t.time());
    CLOG("pointlocation_time: "<<pl_t.time());
    CLOG("rayshooting_time:   "<<rs_t.time());
    CLOG("intersection_time:  "<<it_t.time());
    // warning: the total time showed here could be actually larger
    // that the real time used by this class, since point location
    // and intersection test use the ray shooter and so the same time 
    // could be account to more than one timer
    CLOG("pltotal_time:       "<<
      ct_t.time()+pl_t.time()+rs_t.time()+it_t.time());
  };
};

template <typename SNC_structure>
class SNC_point_locator_by_spatial_subdivision : 
  public SNC_point_locator<SNC_structure>,
  public SNC_decorator<SNC_structure>
{ 
  typedef SNC_decorator<SNC_structure> Base;
  typedef SNC_decorator<SNC_structure> SNC_decorator;
  typedef SNC_point_locator<SNC_structure> SNC_point_locator;
  typedef SNC_point_locator_by_spatial_subdivision<SNC_structure> Self;
  typedef typename SNC_structure::Sphere_map  Sphere_map;
  typedef SM_decorator<Sphere_map>  SM_decorator;
  typedef SNC_intersection<SNC_structure> SNC_intersection;

  typedef typename SNC_structure::Kernel Kernel;

public:
  typedef typename CGAL::SNC_k3_tree_traits<SNC_structure> K3_tree_traits;
  typedef typename CGAL::K3_tree<K3_tree_traits> K3_tree;
  typedef K3_tree SNC_candidate_provider;
  
  #define USING(t) typedef typename SNC_point_locator::t t
  USING(Object_handle);
  USING(Vertex_handle);
  USING(Halfedge_handle);
  USING(Halffacet_handle);
  USING(Halffacet_triangle_handle);
  USING(Vertex_iterator);
  USING(Halfedge_iterator);
  USING(Halffacet_iterator);
  USING(Volume_handle);
  USING(Point_3);
  USING(Segment_3);
  USING(Ray_3);
  USING(Vector_3);
  USING(Triangle_3);
  USING(Aff_transformation_3);
  #undef USING

  typedef typename SNC_candidate_provider::Object_list Object_list;
  typedef typename Object_list::iterator Object_list_iterator;
  typedef typename SNC_candidate_provider::Objects_along_ray Objects_along_ray;
  typedef typename Objects_along_ray::Iterator Objects_along_ray_iterator;

public:
  SNC_point_locator_by_spatial_subdivision() : 
    initialized(false), candidate_provider(0) {}

  virtual void initialize(SNC_structure* W) {
    TIMER(ct_t.start());
    strcpy( version_, "Point Locator by Spatial Subdivision (tm)");
#ifdef CGAL_NEF3_TRIANGULATE_FACETS
    CLOG(version()<<" (with triangulated facets)");
#else
    CLOG(version());
#endif
    CGAL_assertion( W != NULL);
    Base::initialize(W);
    initialized = true;
    Object_list objects;
    Vertex_iterator v;
    Halfedge_iterator e;
    Halffacet_iterator f;
    CGAL_forall_vertices( v, *sncp()) 
      objects.push_back(Object_handle(Vertex_handle(v)));
    CGAL_forall_edges( e, *sncp())
      objects.push_back(Object_handle(Halfedge_handle(e)));
    CGAL_forall_facets( f, *sncp()) {
#ifdef CGAL_NEF3_TRIANGULATE_FACETS
      typedef typename std::list<Triangle_3> Triangles;
      typedef typename Triangles::const_iterator Triangles_iterator;
      typedef Polygon_triangulation_traits_2<Kernel> Triangulation_traits;
      Triangles triangles;
      triangulate_facet<SNC_structure>
	( f, std::back_inserter(triangles), Triangulation_traits());
      for( Triangles_iterator ti = triangles.begin(); 
           ti != triangles.end(); ++ti) {
        Halffacet_triangle_handle th( f, *ti);
        objects.push_back(Object_handle(th));
	CGAL_assertion( assign( th, *(--objects.end())));
	CGAL_assertion( th.get_triangle() == *ti);
      }
#else
      objects.push_back(Object_handle(Halffacet_handle(f)));
#endif // CGAL_NEF3_TRIANGULATE_FACETS
    }
    candidate_provider = new SNC_candidate_provider(objects);
    //_TRACEN(*candidate_provider);
    TIMER(ct_t.stop());
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
    TIMER(ct_t.start());
    CGAL_assertion( initialized);
    bool updated = candidate_provider->update( V, E, F);
    TIMER(ct_t.stop());
    return updated;
  }

  virtual ~SNC_point_locator_by_spatial_subdivision() {
    CGAL_warning( initialized); // required?
    delete candidate_provider;
  }

  virtual Object_handle shoot(const Ray_3& ray) const {
    TIMER(rs_t.start());
    CGAL_assertion( initialized);
    _TRACEN( "shooting: "<<ray);
    Object_handle result;
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    Halffacet_triangle_handle t;
    bool hit = false;
    Point_3 eor; // 'end of ray', the latest ray's hit point
    Objects_along_ray objects = candidate_provider->objects_along_ray(ray);
    Objects_along_ray_iterator objects_iterator = objects.begin();
    while( !hit && objects_iterator != objects.end()) {
      Object_list candidates = *objects_iterator;
      Object_list_iterator o;
      CGAL_for_each( o, candidates) {
        if( assign( v, *o)) {
          _TRACEN("trying vertex on "<<point(v));
          if( ray.source() != point(v) && ray.has_on(point(v))) {
            _TRACEN("the ray intersects the vertex");
            _TRACEN("prev. intersection? "<<hit);
            if( hit) _TRACEN("prev. intersection on "<<eor);
            if( hit && !Segment_3( ray.source(), eor).has_on(point(v)))
              continue;
            eor = point(v);
            result = Object_handle(v);
            hit = true;
            _TRACEN("the vertex becomes the new hit object");
          }
        }
        else if( assign( e, *o)) {
          Point_3 q;
          _TRACEN("trying edge on "<<segment(e));
          if( is.does_intersect_internally( ray, segment(e), q)) {
            _TRACEN("ray intersects edge on "<<q);
            _TRACEN("prev. intersection? "<<hit);
            if( hit) _TRACEN("prev. intersection on "<<eor);
            if( hit && !has_smaller_distance_to_point( ray.source(), q, eor))
              continue;
            _TRACEN("is the intersection point on the current cell? "<<
                    candidate_provider->is_point_on_cell( q, objects_iterator));
            if( !candidate_provider->is_point_on_cell( q, objects_iterator))
                continue;
            eor = q; 
            result = Object_handle(e);
            hit = true;
            _TRACEN("the edge becomes the new hit object");
          }
        }
        else if( assign( f, *o)) {
          Point_3 q;
          _TRACEN("trying facet with on plane "<<plane(f)<<
                  " with point on "<<plane(f).point());
          if( is.does_intersect_internally( ray, f, q) ) {
            _TRACEN("ray intersects facet on "<<q);
            _TRACEN("prev. intersection? "<<hit);
            if( hit) _TRACEN("prev. intersection on "<<eor);
            if( hit && !has_smaller_distance_to_point( ray.source(), q, eor))
	      continue;
            _TRACEN("is the intersection point on the current cell? "<<
                    candidate_provider->is_point_on_cell( q, objects_iterator));
            if( !candidate_provider->is_point_on_cell( q, objects_iterator))
	      continue;
            eor = q;
            result = Object_handle(f);
            hit = true; 
            _TRACEN("the facet becomes the new hit object");
          }
        }
        else if( assign( t, *o)) {
          Point_3 q;
          Triangle_3 tr = t.get_triangle();
          _TRACEN("trying triangle "<<tr);
          if( is.does_intersect( ray, tr, q)) {
            _TRACEN("ray intersect triangle on "<<q);
            _TRACEN("prev. intersection? "<<hit);
            if( hit) _TRACEN("prev. intersection on "<<eor);
            if( hit && !has_smaller_distance_to_point( ray.source(), q, eor))
              continue;
            _TRACEN("is the intersection point on the boundary of the facet? "<<
                    is.does_contain_on_boundary( t, q));
            if( is.does_contain_on_boundary( t, q))
              continue;
            _TRACEN("is the intersection point on the current cell? "<<
                    candidate_provider->is_point_on_cell(q,objects_iterator));
            if( !candidate_provider->is_point_on_cell( q, objects_iterator))
                continue;
            eor = q;
            result = Object_handle(Halffacet_handle(t));
            hit = true; 
            _TRACEN("the facet becomes the new hit object");
          }
        }
        else
          CGAL_assertion_msg( 0, "wrong handle");
      }
      if(!hit)
        ++objects_iterator;
    }
    TIMER(rs_t.stop());
    return result;
  }

  virtual Object_handle locate( const Point_3& p) const {
    TIMER(pl_t.start());
    CGAL_assertion( initialized);
    _TRACEN( "locate "<<p);
    Object_handle result;
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    Halffacet_triangle_handle t;
    Object_list candidates = candidate_provider->objects_around_point(p);
    Object_list_iterator o = candidates.begin();
    bool found = false;
    while( !found && o != candidates.end()) {
      if( assign( v, *o)) {
        if ( p == point(v)) {
          _TRACEN("found on vertex "<<point(v));
          result = Object_handle(v);
          found = true;
        }
      }
      else if( assign( e, *o)) {
        if ( is.does_contain_internally( segment(e), p) ) {
          _TRACEN("found on edge "<<segment(e));
          result = Object_handle(e);
          found = true;
        }
      }
      else if( assign( f, *o)) {
        if ( is.does_contain_internally( f, p) ) {
          _TRACEN("found on facet...");
          result = Object_handle(f);
          found = true;
        }
      }
      else if( assign( t, *o)) {
        Triangle_3 tr = t.get_triangle();
        if( tr.has_on(p)) {
          _TRACEN("found on triangle "<<tr);
	  if(is.does_contain_on_boundary( t, p)) {
	    _TRACEN("but located on the facet's boundary");
	    continue;
	  }
          result = Object_handle(Halffacet_handle(t));
	  found = true;
        }
      }
      o++;
    }
    if( !found) {
      _TRACEN("point not found in 2-skeleton");
      _TRACEN("shooting ray to determine the volume");
      Ray_3 r( p, Vector_3( -1, 0, 0));
      result = Object_handle(determine_volume(r));
    }
    TIMER(pl_t.stop());
    return result;
  }

  virtual void intersect_with_edges( Halfedge_handle e0,
    const typename SNC_point_locator::Intersection_call_back& call_back) const {
    TIMER(it_t.start());
    CGAL_assertion( initialized);
    _TRACEN( "intersecting edge: "<<&*e0<<' '<<segment(e0));
    Segment_3 s(segment(e0));
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    Halffacet_triangle_handle t;
    Object_list_iterator o;
    Object_list objects = candidate_provider->objects_around_segment(s);
    CGAL_for_each( o, objects) {
      if( assign( v, *o)) {
        /* do nothing */
      }
      else if( assign( e, *o)) {
        Point_3 q;
        if( is.does_intersect_internally( s, segment(e), q)) {
          q = normalized(q);
          call_back( e0, Object_handle(e), q);
          _TRACEN("edge intersects edge "<<' '<<&*e<<segment(e)<<" on "<<q);
        }
      }
      else if( assign( f, *o)) {
        /* do nothing */
      }
      else if( assign( t, *o)) {
        /* do nothing */
      }
      else
        CGAL_assertion_msg( 0, "wrong handle");
    }
    TIMER(it_t.stop());
  }

  virtual void intersect_with_facets( Halfedge_handle e0, 
    const typename SNC_point_locator::Intersection_call_back& call_back) const {
    TIMER(it_t.start());
    CGAL_assertion( initialized);
    _TRACEN( "intersecting edge: "<<segment(e0));
#ifdef CGAL_NEF3_TRIANGULATE_FACETS
    Unique_hash_map< Halffacet_triangle_handle, bool> f_mark(false);
#endif // CGAL_NEF3_TRIANGULATE_FACETS
    Segment_3 s(segment(e0));
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    Halffacet_triangle_handle t;
    Object_list_iterator o;
    Object_list objects = candidate_provider->objects_around_segment(s);
    CGAL_for_each( o, objects) {
      if( assign( v, *o)) {
        /* do nothing */
      }
      else if( assign( e, *o)) {
        /* do nothing */
      }
      else if( assign( f, *o)) {
        Point_3 q;
        if( is.does_intersect_internally( s, f, q) ) {
          q = normalized(q);
          call_back( e0, Object_handle(f), q);
          _TRACEN("edge intersects facet on plane "<<plane(f)<<" on "<<q);
        }
      }
      else if( assign( t, *o)) {
        Point_3 q;
        Triangle_3 tr = t.get_triangle();
#ifdef CGAL_NEF3_TRIANGULATE_FACETS
	if( f_mark[t])
	  continue;
#endif // CGAL_NEF3_TRIANGULATE_FACETS
	_TRACEN("trying with triangle "<<tr);
        if( is.does_intersect( s, tr, q) &&
            !is.does_contain_on_boundary( t, q)) {
          q = normalized(q);
          call_back( e0, Object_handle(Halffacet_handle(t)), q);
          _TRACEN("edge intersects facet triangle on plane "<<plane(t)<<" on "<<q);
#ifdef CGAL_NEF3_TRIANGULATE_FACETS
	  f_mark[t] = true;
#endif // CGAL_NEF3_TRIANGULATE_FACETS
        }
      }
      else
        CGAL_assertion_msg( 0, "wrong handle");
    }
    TIMER(it_t.stop());
  }

private:
  Volume_handle determine_volume( const Ray_3& ray) const {
    Halffacet_handle f_below;
    Object_handle o = shoot(ray);
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    if( assign( v, o)) {
      _TRACEN("vertex hit, obtaining volume...");
      f_below = get_visible_facet( v, ray);
      if( f_below != Halffacet_handle())
        return volume(f_below);
      SM_decorator SD(&*v); // now, the vertex has no incident facets
      CGAL_assertion( SD.number_of_sfaces() == 1);
      return volume(SD.sfaces_begin());
    }
    else if( assign( e, o)) {
      _TRACEN("edge hit, obtaining volume...");
      f_below = get_visible_facet( e, ray);
      if( f_below != Halffacet_handle())
        return volume(f_below);
      CGAL_assertion_code(SM_decorator SD(&*source(e))); // now, the edge has no incident facets
      CGAL_assertion(SD.is_isolated(e));
      return volume(sface(e));
    }
    else if( assign( f, o)) {
      _TRACEN("facet hit, obtaining volume...");
      f_below = get_visible_facet(f, ray);
      CGAL_assertion( f_below != Halffacet_handle());
      return volume(f_below);
    }
    return Base(*this).volumes_begin(); // TODO: Comment this hack!
  }

private:
  bool initialized;
  SNC_candidate_provider* candidate_provider;
  SNC_intersection is;

  std::list<Halffacet_triangle_handle> triangulation;
};

template <typename SNC_structure>
class SNC_point_locator_naive : 
  public SNC_ray_shooter<SNC_structure>, 
  public SNC_point_locator<SNC_structure>
{
  typedef SNC_ray_shooter<SNC_structure> Base;
  typedef SNC_point_locator_naive<SNC_structure> Self;
  typedef SNC_point_locator<SNC_structure> SNC_point_locator;
  typedef SNC_intersection<SNC_structure> SNC_intersection;

public:
  #define USING(t) typedef typename SNC_point_locator::t t
  USING(Object_handle);
  USING(Vertex_handle);
  USING(Halfedge_handle);
  USING(Halffacet_handle);
  USING(Volume_handle);
  USING(Vertex_iterator);
  USING(Halfedge_iterator);
  USING(Halffacet_iterator);
  USING(Point_3);
  USING(Segment_3);
  USING(Ray_3);
  USING(Triangle_3);
  #undef USING

public:
  SNC_point_locator_naive() : initialized(false) {}
  virtual void initialize(SNC_structure* W) { 
    TIMER(ct_t.start());
    strcpy( version_, "Naive Point Locator (tm)");
    CLOG(version());
    CGAL_assertion( W != NULL);
    Base::initialize(W); 
    initialized = true;
    TIMER(ct_t.stop());
  }

  virtual Self* clone() const { 
    return new Self; 
  }

  virtual bool update( Unique_hash_map<Vertex_handle, bool>& V, 
                       Unique_hash_map<Halfedge_handle, bool>& E, 
                       Unique_hash_map<Halffacet_handle, bool>& F) {
    TIMER(ct_t.start());
    CGAL_assertion( initialized);
    TIMER(ct_t.stop());
    return false;
  }

  virtual ~SNC_point_locator_naive() {}

  virtual Object_handle locate(const Point_3& p) const {
    TIMER(pl_t.start());
    CGAL_assertion( initialized);
    TIMER(pl_t.stop());
    return Base::locate(p);
  }

  virtual Object_handle shoot(const Ray_3& r) const {
    TIMER(rs_t.start());
    CGAL_assertion( initialized);
    TIMER(rs_t.stop());
    return Base::shoot(r);
  }

  virtual void intersect_with_edges( Halfedge_handle e0, 
    const typename SNC_point_locator::Intersection_call_back& call_back) const {
    TIMER(it_t.start());
    CGAL_assertion( initialized);
    TRACEN( "intersecting edge: "<<segment(e0));
    SNC_intersection is(*sncp());
    Segment_3 s(segment(e0));
    Halfedge_iterator e;
    CGAL_forall_edges( e, *sncp()) {
      Point_3 q;
      if( is.does_intersect_internally( s, segment(e), q)) {
        q = normalized(q);
        TRACEN("edge intersects edge "<<segment(e)<<" on "<<q);
        call_back( e0, Object_handle(e), q);
      }
    }
    TIMER(it_t.stop());
  }

  virtual void intersect_with_facets( Halfedge_handle e0, 
    const typename SNC_point_locator::Intersection_call_back& call_back) const {
    TIMER(it_t.start());
    CGAL_assertion( initialized);
    TRACEN( "intersecting edge: "<<segment(e0));
    SNC_intersection is(*sncp());
    Segment_3 s(segment(e0));
    Halffacet_iterator f;
    CGAL_forall_facets( f, *sncp()) {
      Point_3 q;
      if( is.does_intersect_internally( s, f, q) ) {
        q = normalized(q);
        TRACEN("edge intersects facet on plane "<<plane(f)<<" on "<<q);
        call_back( e0, Object_handle(f), q);
      }
    }
    TIMER(it_t.stop());
  }

private:
  bool initialized;
};


CGAL_END_NAMESPACE


#endif // SNC_POINT_LOCATOR_H

