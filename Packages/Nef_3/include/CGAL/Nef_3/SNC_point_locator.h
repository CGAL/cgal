#line 7 "point_locator.nw"
#ifndef SNC_POINT_LOCATOR_H
#define SNC_POINT_LOCATOR_H

#include <CGAL/basic.h>
#include <CGAL/Nef_3/SNC_intersection.h>
#include <CGAL/Nef_3/SNC_ray_shooter.h>
#include <CGAL/Nef_3/SNC_k3_tree_traits.h>
#include <CGAL/Nef_3/K3_tree.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Timer.h>

// #include <CGAL/Polygon_triangulation_traits_2.h>
// #include <CGAL/Nef_3/triangulate_nef3_facet.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 509
#include <CGAL/Nef_2/debug.h>

#undef _TRACEN
#define _TRACEN(msg) TRACEN( "SNC_point_locator: " << msg);

// TODO: find out the proper CGAL replacement for this macro and remove it
#define CGAL_for_each( i, C) for( i = C.begin(); i != C.end(); ++i)

// #define TIMER(instruction) instruction
#define TIMER(instruction)

// #define CLOG(t) std::clog <<" "<<t<<std::endl; std::clog.flush()
#define CLOG(t)

CGAL_BEGIN_NAMESPACE

template <typename SNC_decorator>
class SNC_point_locator
{
 public:
  class Intersection_call_back;
  typedef SNC_point_locator<SNC_decorator> Self;
  typedef typename SNC_decorator::Decorator_traits Decorator_traits;
  typedef typename SNC_decorator::SNC_structure SNC_structure;
protected:
  char version_[64];
  // time for construction, point location, ray shooting and intersection test
  mutable Timer ct_t, pl_t, rs_t, it_t; 

public: 
  typedef typename SNC_structure::Object_handle Object_handle;
  typedef typename SNC_structure::Halfedge_const_handle Halfedge_const_handle;
  typedef typename SNC_structure::Halffacet_const_handle Halffacet_const_handle;
  typedef typename SNC_structure::Halffacet_triangle_handle 
                                  Halffacet_triangle_handle;
  typedef typename SNC_structure::Point_3 Point_3;
  typedef typename SNC_structure::Segment_3 Segment_3;
  typedef typename SNC_structure::Ray_3 Ray_3;
  typedef typename SNC_structure::Vector_3 Vector_3;
  typedef typename SNC_structure::Triangle_3 Triangle_3;
  typedef typename SNC_structure::Aff_transformation_3 
                                  Aff_transformation_3;

  typedef typename Decorator_traits::Vertex_handle Vertex_handle;
  typedef typename Decorator_traits::Halfedge_handle Halfedge_handle;
  typedef typename Decorator_traits::Halffacet_handle Halffacet_handle;
  typedef typename Decorator_traits::Volume_handle Volume_handle;
  typedef typename Decorator_traits::Vertex_iterator Vertex_iterator;
  typedef typename Decorator_traits::Halfedge_iterator Halfedge_iterator;
  typedef typename Decorator_traits::Halffacet_iterator Halffacet_iterator;

  const char* version() { return version_; }

  virtual Object_handle locate(const Point_3& p) const = 0;

  virtual Object_handle shoot(const Ray_3& s) const = 0;

  virtual void intersect_with_edges( Halfedge_const_handle edge,
                                     const Intersection_call_back& call_back) 
    const = 0;

  virtual void intersect_with_facets( Halfedge_const_handle edge,
                                      const Intersection_call_back& call_back)
    const = 0;

  virtual void intersect_with_edges_and_facets( Halfedge_const_handle edge,
	const Intersection_call_back& call_back) const = 0;

  class Intersection_call_back 
  {
  public:
    virtual void operator()( Halfedge_const_handle edge, Object_handle object, 
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
  typedef typename SNC_structure::Halfedge_const_handle Halfedge_const_handle;
  typedef typename SNC_structure::Halffacet_const_handle Halffacet_const_handle; 
  typedef typename SNC_structure::Halffacet_triangle_handle 
                                  Halffacet_triangle_handle;	
  typedef typename SNC_structure::Point_3 Point_3;
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

  typedef typename SNC_candidate_provider::Object_list Object_list;
  typedef typename Object_list::iterator Object_list_iterator;
  typedef typename SNC_candidate_provider::Objects_along_ray Objects_along_ray;
  typedef typename Objects_along_ray::Iterator Objects_along_ray_iterator;

public:
  SNC_point_locator_by_spatial_subdivision() : 
    initialized(false), candidate_provider(0) {}

  virtual void initialize(SNC_structure* W) {
    TIMER(ct_t.start());
    strcpy( this->version_, "Point Locator by Spatial Subdivision (tm)");
#ifdef CGAL_NEF3_TRIANGULATE_FACETS
    CLOG(version()<<" (with triangulated facets)");
#else
    CLOG(version());
#endif
    CGAL_assertion( W != NULL);
//    (Base) *this = SNC_decorator(*W);
	set_snc(*W);
    initialized = true;
    Object_list objects;
    Vertex_iterator v;
    Halfedge_iterator e;
    Halffacet_iterator f;
    CGAL_forall_vertices( v, *this->sncp())
      objects.push_back(Object_handle(Vertex_handle(v)));
    typename Object_list::difference_type v_end = objects.size();
    CGAL_forall_edges( e, *this->sncp())
      objects.push_back(Object_handle(Halfedge_handle(e)));
    CGAL_forall_facets( f, *this->sncp()) {
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
	CGAL_assertion( CGAL::assign( th, *(--objects.end())));
	CGAL_assertion( th.get_triangle() == *ti);
      }
#else
      objects.push_back(Object_handle(Halffacet_handle(f)));
#endif // CGAL_NEF3_TRIANGULATE_FACETS
    }
    Object_list_iterator oli=objects.begin()+v_end;
    candidate_provider = new SNC_candidate_provider(objects,oli);
//    TRACEN(*candidate_provider);
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
        if( CGAL::assign( v, *o)) {
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
        else if( CGAL::assign( e, *o)) {
          Point_3 q;
          _TRACEN("trying edge on "<< Segment_3(e->source()->point(),e->twin()->source()->point()));
          if( is.does_intersect_internally( ray, Segment_3(e->source()->point(),
                                                           e->twin()->source()->point()), q)) {
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
        else if( CGAL::assign( f, *o)) {
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
        else if( CGAL::assign( t, *o)) {
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
    if(Infi_box::extended_kernel()) {
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
      if( CGAL::assign( v, *o)) {
        if ( p == point(v)) {
          _TRACEN("found on vertex "<<point(v))          
          result = Object_handle(v);
          found = true;
        }
      }
      else if( CGAL::assign( e, *o)) {
        if ( is.does_contain_internally(Segment_3(e->source()->point(),e->twin()->source()->point()), p) ) {
          _TRACEN("found on edge "<<Segment_3(e->source()->point(),e->twin()->source()->point()));
          result = Object_handle(e);
          found = true;
        }
      }
      else if( CGAL::assign( f, *o)) {
        if ( is.does_contain_internally( f, p) ) {
          _TRACEN("found on facet...");
          result = Object_handle(f);
          found = true;
        }
      }
      else if( CGAL::assign( t, *o)) {
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
    }    TIMER(pl_t.start());
    TIMER(pl_t.stop());
    return result;
  } else {   // standard kernel
    CGAL_assertion( initialized);
    _TRACEN( "locate "<<p);
    typename SNC_structure::FT min_distance;
    typename SNC_structure::FT tmp_distance;
    Object_handle result;
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    Halffacet_triangle_handle t;
    Object_list candidates = candidate_provider->objects_around_point(p);
    Object_list_iterator o = candidates.begin();

    if(candidates.empty())
      return Base(*this).volumes_begin();

    CGAL::assign(v,*o);
    CGAL_assertion(CGAL::assign(v,*o));
    if(p==point(v))
      return Object_handle(v);

    min_distance = CGAL::squared_distance(point(v),p);
    result = Object_handle(v);
    ++o;
    while(o!=candidates.end() && CGAL::assign(v,*o)) {
      if ( p == point(v)) {
        _TRACEN("found on vertex "<<point(v));
        return Object_handle(v);
      }
      tmp_distance = CGAL::squared_distance(v->point(),p);
      if(tmp_distance < min_distance) {
	result = Object_handle(v);
        min_distance = tmp_distance;
      }
      ++o;
    }     

    CGAL::assign(v, result);
    Segment_3 s(p,point(v));
    Point_3 ip;

    Object_list_iterator ox(o);
    for(;o!=candidates.end();++o) {
      if( CGAL::assign( e, *o)) {
	Segment_3 ss(e->source()->point(),e->twin()->source()->point());
	TRACEN("test edge " << e->source()->point() << "->" << e->twin()->source()->point());
        if(is.does_contain_internally(ss, p) ) {
          _TRACEN("found on edge "<< ss);
          return Object_handle(e);
        }
	if(is.does_intersect_internally(s, ss, ip)) {
	  s = Segment_3(p, normalized(ip));
	  result = Object_handle(e);
        }
      } else if( CGAL::assign( f, *o)) {
	TRACEN("test facet " << f->plane());
        if (is.does_contain_internally(f,p) ) {
          _TRACEN("found on facet...");
          return Object_handle(f);
        }
        if( is.does_intersect_internally(s,f,ip)) {	
          s = Segment_3(p, normalized(ip));
	  result = Object_handle(f);
        }
      } else CGAL_assertion_msg(false, "wrong handle type");
    }

    if( CGAL::assign( v, result)) {
      _TRACEN("vertex hit, obtaining volume...");
      SM_point_locator L(&*v);
      Object_handle so = L.locate(s.source()-s.target());
      SFace_handle sf;
      if(CGAL::assign(sf,so))
        return sf->volume();

      CGAL_assertion_msg(false, "wrong handle type");
/*
      SHalfedge_handle se;
      CGAL_assertion(CGAL::assign(se,so));
      TRACEN("intersect segment " << s << " with edges");
      for(;ox!=candidates.end();++ox) {
	if(!CGAL::assign(e,*ox)) continue;
	TRACEN("test edge " << e->source()->point() << "->" << e->twin()->source()->point());
	if(is.does_intersect_internally(s,Segment_3(e->source()->point(),e->twin()->source()->point()),ip)) {
	  s = Segment_3(p, normalized(ip));
	  result = Object_handle(e);
        }
      }
      CGAL_assertion(CGAL::assign(e,result));
      CGAL::assign(e,result);
      f = get_visible_facet(e, Ray_3(p, s.target()));	
      if( f != Halffacet_handle())
	return f->incident_volume();
      SM_decorator SD(&*v); // now, the vertex has no incident facets
      CGAL_assertion( SD.number_of_sfaces() == 1);
      return volume(SD.sfaces_begin());      
*/
    } else if( CGAL::assign( f, result)) {
      _TRACEN("facet hit, obtaining volume...");
      if(f->plane().oriented_side(p) == ON_NEGATIVE_SIDE)
	f = f->twin();
      return Object_handle(f->incident_volume());
    } else if( CGAL::assign(e, result)) {
      SM_decorator SD(&*source(e));
      if( SD.is_isolated(e))
        return Object_handle(e->incident_sface()->volume());	
      return get_visible_facet(e,Ray_3(s.source(),s.to_vector()))->incident_volume();
    }
    CGAL_assertion_msg(false, "wrong handle type");
    return Object_handle();
  }
  }

  virtual void intersect_with_edges_and_facets( Halfedge_const_handle e0,
	const typename SNC_point_locator::Intersection_call_back& call_back) const {

    TIMER(it_t.start());
    CGAL_assertion( initialized);
    _TRACEN( "intersecting edge: "<<&*e0<<' '<<Segment_3(e0->source()->point(),
                                                         e0->twin()->source()->point()));

#ifdef CGAL_NEF3_TRIANGULATE_FACETS
    Unique_hash_map< Halffacet_triangle_handle, bool> f_mark(false);
#endif // CGAL_NEF3_TRIANGULATE_FACETS

    Segment_3 s(Segment_3(e0->source()->point(),e0->twin()->source()->point()));
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    Halffacet_triangle_handle t;
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
          call_back( e0, Object_handle(Halfedge_const_handle(e)), q);
          _TRACEN("edge intersects edge "<<' '<<&*e<< Segment_3(e->source()->point(),
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
          call_back( e0, Object_handle(Halffacet_const_handle(f)), q);
          _TRACEN("edge intersects facet on plane "<<plane(f)<<" on "<<q);
        }
      }
      else if( CGAL::assign( t, *o)) {
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

  virtual void intersect_with_edges( Halfedge_const_handle e0,
    const typename SNC_point_locator::Intersection_call_back& call_back) const {
    TIMER(it_t.start());
    CGAL_assertion( initialized);
    _TRACEN( "intersecting edge: "<<&*e0<<' '<<Segment_3(e0->source()->point(),
                                                         e0->twin()->source()->point()));
    Segment_3 s(Segment_3(e0->source()->point(),e0->twin()->source()->point()));
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    Halffacet_triangle_handle t;
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
          call_back( e0, Object_handle(Halfedge_const_handle(e)), q);
          _TRACEN("edge intersects edge "<<' '<<&*e<< Segment_3(e->source()->point(),
                                                                e->twin()->source()->point())<<" on "<<q);
        }
      }
      else if( CGAL::assign( f, *o)) {
        /* do nothing */
      }
      else if( CGAL::assign( t, *o)) {
        /* do nothing */
      }
      else
        CGAL_assertion_msg( 0, "wrong handle");
    }
    TIMER(it_t.stop());
  }

  virtual void intersect_with_facets( Halfedge_const_handle e0, 
    const typename SNC_point_locator::Intersection_call_back& call_back) const {
    TIMER(it_t.start());
    CGAL_assertion( initialized);
    _TRACEN( "intersecting edge: "<< Segment_3(e0->source()->point(),
                                               e0->twin()->source()->point()));
#ifdef CGAL_NEF3_TRIANGULATE_FACETS
    Unique_hash_map< Halffacet_triangle_handle, bool> f_mark(false);
#endif // CGAL_NEF3_TRIANGULATE_FACETS
    Segment_3 s(Segment_3(e0->source()->point(),e0->twin()->source()->point()));
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    Halffacet_triangle_handle t;
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
          call_back( e0, Object_handle(Halffacet_const_handle(f)), q);
          _TRACEN("edge intersects facet on plane "<<plane(f)<<" on "<<q);
        }
      }
      else if( CGAL::assign( t, *o)) {
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
    if( CGAL::assign( v, o)) {
      _TRACEN("vertex hit, obtaining volume...");
      f_below = get_visible_facet( v, ray);
      if( f_below != Halffacet_handle())
        return f_below->incident_volume();
      SM_decorator SD(&*v); // now, the vertex has no incident facets
      CGAL_assertion( SD.number_of_sfaces() == 1);
      return volume(SD.sfaces_begin());
    }
    else if( CGAL::assign( e, o)) {
      _TRACEN("edge hit, obtaining volume...");
      f_below = get_visible_facet( e, ray);
      if( f_below != Halffacet_handle())
        return f_below->incident_volume();
      CGAL_assertion_code(SM_decorator SD(&*source(e))); // now, the edge has no incident facets
      CGAL_assertion(SD.is_isolated(e));
      return volume(sface(e));
    }
    else if( CGAL::assign( f, o)) {
      _TRACEN("facet hit, obtaining volume...");
      f_below = get_visible_facet(f, ray);
      CGAL_assertion( f_below != Halffacet_handle());
      return f_below->incident_volume();
    }
    return Base(*this).volumes_begin(); // TODO: Comment this hack!
  }

private:
  bool initialized;
  SNC_candidate_provider* candidate_provider;
  SNC_intersection is;
 
 std::list<Halffacet_triangle_handle> triangulation;
};

/*
template <typename SNC_decorator>
class SNC_point_locator_naive : 
  public SNC_ray_shooter<SNC_structure>, 
  public SNC_point_locator<SNC_structure>
{
  typedef SNC_ray_shooter<SNC_structure> Base;
  typedef SNC_point_locator_naive<SNC_structure> Self;
  typedef SNC_point_locator<SNC_structure> SNC_point_locator;
  typedef SNC_intersection<SNC_structure> SNC_intersection;

public:
  typedef typename SNC_decorator::Object_handle Object_handle;
  typedef typename SNC_decorator::Halfedge_const_handle Halfedge_const_handle;
  typedef typename SNC_decorator::Halffacet_triangle_handle 
                                  Halffacet_triangle_handle;
  typedef typename SNC_decorator::Point_3 Point_3;
  typedef typename SNC_decorator::Segment_3 Segment_3;
  typedef typename SNC_decorator::Ray_3 Ray_3;
  typedef typename SNC_decorator::Triangle_3 Triangle_3;


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

  virtual void intersect_with_edges( Halfedge_const_handle e0, 
    const typename SNC_point_locator::Intersection_call_back& call_back) const {
    TIMER(it_t.start());
    CGAL_assertion( initialized);
    TRACEN( "intersecting edge: "<< Segment_3(e0->source()->point(),
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
        TRACEN("edge intersects edge "<< Segment_3(e->source()->point(),
                                                   e->twin()->source()->point()) <<" on "<<q);
        call_back( e0, Object_handle(e), q);
      }
    }
    TIMER(it_t.stop());
  }

  virtual void intersect_with_facets( Halfedge_const_handle e0, 
    const typename SNC_point_locator::Intersection_call_back& call_back) const {
    TIMER(it_t.start());
    CGAL_assertion( initialized);
    TRACEN( "intersecting edge: "<< Segment_3(e0->source()->point(),
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
        TRACEN("edge intersects facet on plane "<<plane(f)<<" on "<<q);
        call_back( e0, Object_handle(f), q);
      }
    }
    TIMER(it_t.stop());
  }

private:
  bool initialized;
};

*/

CGAL_END_NAMESPACE
#endif // SNC_POINT_LOCATOR_H

