#line 7 "point_locator.nw"
#ifndef SNC_POINT_LOCATOR_H
#define SNC_POINT_LOCATOR_H

#include <CGAL/basic.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_SM_point_locator.h>
#include <CGAL/Nef_3/SNC_intersection.h>
#include <CGAL/Nef_3/SNC_ray_shooter.h>
#include <CGAL/Nef_3/SNC_k3_tree_traits.h>
#include <CGAL/K3_tree.h>
#include <CGAL/Unique_hash_map.h>

#undef _DEBUG
#define _DEBUG 509
#include <CGAL/Nef_3/debug.h>

#undef _TRACEN
#define _TRACEN(msg) TRACEN( "SNC_point_locator: " << msg);

#define CGAL_for_each( i, C) for( i = C.begin(); i != C.end(); ++i)
// TODO: find out the CGAL replacement for this macro and remove it

CGAL_BEGIN_NAMESPACE

template <typename SNC_structure>
class SNC_point_locator
{
  class Intersection_call_back;
  typedef SNC_point_locator<SNC_structure> Self;
  typedef typename SNC_structure::Object_handle Object_handle;
  typedef typename SNC_structure::Vertex_handle Vertex_handle;
  typedef typename SNC_structure::Halfedge_handle Halfedge_handle;
  typedef typename SNC_structure::Halffacet_handle Halffacet_handle;
  typedef typename SNC_structure::Point_3 Point_3;
  typedef typename SNC_structure::Segment_3 Segment_3;
  typedef typename SNC_structure::Ray_3 Ray_3;
public:
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
  //virtual Self* copy() const = 0;
  virtual bool update( Unique_hash_map<Vertex_handle, bool>& V, 
                       Unique_hash_map<Halfedge_handle, bool>& E, 
                       Unique_hash_map<Halffacet_handle, bool>& F) = 0;
  virtual ~SNC_point_locator() {};
};

template <typename SNC_structure>
class SNC_point_locator_by_spatial_subdivision : 
  public SNC_point_locator<SNC_structure>,
  public SNC_decorator<SNC_structure>
{ 
  typedef SNC_decorator<SNC_structure> Base;
  typedef SNC_point_locator<SNC_structure> SNC_point_locator;
  typedef SNC_decorator<SNC_structure> SNC_decorator;
  typedef SNC_point_locator_by_spatial_subdivision<SNC_structure> Self;
  typedef SNC_SM_decorator<SNC_structure> SM_decorator;
  typedef SNC_intersection<SNC_structure> SNC_intersection;

  typedef typename CGAL::SNC_k3_tree_traits<SNC_structure> K3_tree_traits;
  typedef typename CGAL::K3_tree<K3_tree_traits> K3_tree;
  typedef K3_tree SNC_candidate_provider;

  typedef typename SNC_candidate_provider::Object_list Object_list;
  typedef typename Object_list::iterator Object_list_iterator;
  typedef typename SNC_candidate_provider::Objects_along_ray Objects_along_ray;
  typedef typename Objects_along_ray::Iterator Objects_along_ray_iterator;

  #define USING(t) typedef typename SNC_structure::t t
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
  USING(Direction_3);
  #undef USING

  bool initialized;

public:
  SNC_point_locator_by_spatial_subdivision() : initialized(false) {}
  virtual void initialize(SNC_structure* W) {
    CGAL_assertion( W != NULL);
    Base::initialize(W);
    initialized = true;
    Object_list objects;
    Vertex_iterator v;
    Halfedge_iterator e;
    Halffacet_iterator f;
    CGAL_nef3_forall_vertices( v, *sncp()) 
      objects.push_back(Object_handle(Vertex_handle(v)));
    CGAL_nef3_forall_edges( e, *sncp())
      objects.push_back(Object_handle(Halfedge_handle(e)));
    CGAL_nef3_forall_facets( f, *sncp())
      objects.push_back(Object_handle(Halffacet_handle(f)));
    candidate_provider = new SNC_candidate_provider(objects);
    _TRACEN(*candidate_provider);
  }

  virtual Self* clone() const { 
    return new Self; 
  }
  /*
  virtual Self* copy() const { 
    CGAL_assertion( initialized);
    Self *pl_copy = new Self;
    pl_copy->initialize(sncp());
    return pl_copy; 
  }
  */
  virtual bool update( Unique_hash_map<Vertex_handle, bool>& V, 
                       Unique_hash_map<Halfedge_handle, bool>& E, 
                       Unique_hash_map<Halffacet_handle, bool>& F) {
    CGAL_assertion( initialized);
    _TRACEN( *candidate_provider);
    bool updated = candidate_provider->update( V, E, F);
    _TRACEN( *candidate_provider);
    return updated;
  }
  virtual ~SNC_point_locator_by_spatial_subdivision() {
    CGAL_assertion( initialized);
    delete candidate_provider;
  }

  virtual Object_handle shoot(const Ray_3& ray) const {
    CGAL_assertion( initialized);
    _TRACEN( "shooting: "<<ray);
    Object_handle result;
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    bool hit = false;
    Point_3 end_of_seg;
    Objects_along_ray objects = candidate_provider->objects_along_ray(ray);
    Objects_along_ray_iterator objects_iterator = objects.begin();
    while( !hit && objects_iterator != objects.end()) {
      Object_list candidates = *objects_iterator;
      Object_list_iterator o;
      CGAL_for_each( o, candidates) {
        if( assign( v, *o)) {
	  _TRACEN("trying vertex on "<<point(v));
          if( ray.source() != point(v) && ray.has_on(point(v))) {
	    _TRACEN("the ray intersects the point");
	    _TRACEN("actual intersection? "<<hit<<" on "<<end_of_seg);
            if( hit && 
                !Segment_3( ray.source(), end_of_seg).has_on(point(v)))
              continue;
            end_of_seg = point(v);
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
            _TRACEN("actual intersection? "<<hit<<" on "<<end_of_seg);
	    _TRACEN("is the intersection point on the current cell? "<<candidate_provider->is_point_on_cell( q, objects_iterator));
	    if( !candidate_provider->is_point_on_cell( q, objects_iterator)) // TODO: hide the node
		continue;
            if( !hit || 
		has_smaller_distance_to_point( ray.source(), q, end_of_seg)) {
	      end_of_seg = q; 
	      result = Object_handle(e);
	      hit = true;
	      _TRACEN("the edge becomes the new hit object");
	    }
          }
        }
        else if( assign( f, *o)) {
          Point_3 q;
	  _TRACEN("trying facet with on plane "<<plane(f)<<" with point on "<<plane(f).point());
          if( is.does_intersect_internally( ray, f, q) ) {
	    _TRACEN("ray intersects facet on "<<q);
            _TRACEN("actual intersection? "<<hit<<" on "<<end_of_seg);
	    _TRACEN("is the intersection point on the current cell? "<<candidate_provider->is_point_on_cell( q, objects_iterator));
	    if( !candidate_provider->is_point_on_cell( q, objects_iterator))
		continue;
            if( !hit || 
		has_smaller_distance_to_point( ray.source(), q, end_of_seg)) {
	      end_of_seg = q;
	      result = Object_handle(f);
	      hit = true; 
	      _TRACEN("the facet becomes the new hit object");
	    }
          }
        }
        else
          CGAL_nef3_assertion_msg( 0, "wrong handle");
      }
      if(!hit)
	++objects_iterator;
    }
    return result;
  }

  virtual Object_handle locate( const Point_3& p) const {
    CGAL_assertion( initialized);
    _TRACEN( "locate "<<p);
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    Object_list candidates = candidate_provider->objects_around_point(p);
    Object_list_iterator o;
    CGAL_for_each( o, candidates) {
      if( assign( v, *o)) {
        if ( p == point(v)) {
          _TRACEN("found on vertex"<<point(v));
          return Object_handle(v);
        }
      }
      else if( assign( e, *o)) {
        if ( is.does_contain_internally( segment(e), p) ) {
          _TRACEN("found on edge"<<segment(e));
          return Object_handle(e);
        }
      }
      else if( assign( f, *o)) {
        if ( is.does_contain_internally( f, p) ) {
          _TRACEN("found on facet...");
          return Object_handle(f);
        }
      }
    }
    _TRACEN("point not found in 2-skeleton");

    // TODO: here, the vertex to choose should be contained in the same cell 
    // where the query point is located, or, if there is not vertex on the 
    // cell, then the volume containing the cell should be available in the 
    // tree structure and so immediately reported
    _TRACEN("shooting ray to determine the volume");
    Ray_3 r( p, Direction_3( -1, 0, 0));
    return Object_handle(determine_volume(r));
  }

  virtual void intersect_with_edges( Halfedge_handle e0,
    const typename SNC_point_locator::Intersection_call_back& call_back) const {
    CGAL_assertion( initialized);
    _TRACEN( "intersecting edge: "<<&*e0<<' '<<segment(e0));
    Segment_3 s(segment(e0));
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
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
      else
        CGAL_nef3_assertion_msg( 0, "wrong handle");
    }
  }

  void intersect_with_facets( Halfedge_handle e0, 
    const typename SNC_point_locator::Intersection_call_back& call_back) const {
    CGAL_assertion( initialized);
    _TRACEN( "intersecting edge: "<<segment(e0));
    Segment_3 s(segment(e0));
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
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
      else
        CGAL_nef3_assertion_msg( 0, "wrong handle");
    }
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
      SM_decorator SD(v); // now, the vertex has no incident facets
      CGAL_nef3_assertion( SD.number_of_sfaces() == 1);
      return volume(SD.sfaces_begin());
    }
    else if( assign( e, o)) {
      _TRACEN("edge hit, obtaining volume...");
      f_below = get_visible_facet( e, ray);
      if( f_below != Halffacet_handle())
        return volume(f_below);
      SM_decorator SD(source(e)); // now, the edge has no incident facets
      CGAL_nef3_assertion(SD.is_isolated(e));
      return volume(sface(e));
    }
    else if( assign( f, o)) {
      _TRACEN("facet hit, obtaining volume...");
      f_below = get_visible_facet(f, ray);
      CGAL_nef3_assertion( f_below != Halffacet_handle());
      return volume(f_below);
    }
    return Base(*this).volumes_begin(); // TODO: Comment this hack!
  }

private:
  SNC_candidate_provider* candidate_provider;
  SNC_intersection is;
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
  typedef typename SNC_structure::Object_handle Object_handle;
  typedef typename SNC_structure::Point_3 Point_3;
  typedef typename SNC_structure::Segment_3 Segment_3;
  typedef typename SNC_structure::Ray_3 Ray_3;
  typedef typename SNC_structure::Vertex_handle Vertex_handle;
  typedef typename SNC_structure::Halfedge_handle Halfedge_handle;
  typedef typename SNC_structure::Halffacet_handle Halffacet_handle;
  typedef typename SNC_structure::Vertex_iterator Vertex_iterator;
  typedef typename SNC_structure::Halfedge_iterator Halfedge_iterator;
  typedef typename SNC_structure::Halffacet_iterator Halffacet_iterator;

  bool initialized;

public:
  SNC_point_locator_naive() : initialized(false) {}
  virtual void initialize(SNC_structure* W) { 
    CGAL_assertion( W != NULL);
    Base::initialize(W); 
    initialized = true;
  }

  virtual Self* clone() const { 
    return new Self; 
  }
  /*
  virtual Self* copy() const { 
    CGAL_assertion( initialized);
    Self *pl_copy = new Self;
    pl_copy->initialize(sncp());
    return pl_copy; 
  }
  */
  virtual bool update( Unique_hash_map<Vertex_handle, bool>& V, 
                       Unique_hash_map<Halfedge_handle, bool>& E, 
                       Unique_hash_map<Halffacet_handle, bool>& F) {
    CGAL_assertion( initialized);
    return false;
  }

  virtual ~SNC_point_locator_naive() {}

  Object_handle locate(const Point_3& p) const {
    CGAL_assertion( initialized);
    return Base::locate(p);
  }
  Object_handle shoot(const Ray_3& r) const {
    CGAL_assertion( initialized);
    return Base::shoot(r);
  }

  void intersect_with_edges( Halfedge_handle e0, 
    const typename SNC_point_locator::Intersection_call_back& call_back) const {
    CGAL_assertion( initialized);
    TRACEN( "intersecting edge: "<<segment(e0));
    SNC_intersection is(*sncp());
    Segment_3 s(segment(e0));
    Halfedge_iterator e;
    CGAL_nef3_forall_edges( e, *sncp()) {
      Point_3 q;
      if( is.does_intersect_internally( s, segment(e), q)) {
        q = normalized(q);
        TRACEN("edge intersects edge "<<segment(e)<<" on "<<q);
        call_back( e0, Object_handle(e), q);
      }
    }
  }

  void intersect_with_facets( Halfedge_handle e0, 
    const typename SNC_point_locator::Intersection_call_back& call_back) const {
    CGAL_assertion( initialized);
    TRACEN( "intersecting edge: "<<segment(e0));
    SNC_intersection is(*sncp());
    Segment_3 s(segment(e0));
    Halffacet_iterator f;
    CGAL_nef3_forall_facets( f, *sncp()) {
      Point_3 q;
      if( is.does_intersect_internally( s, f, q) ) {
        q = normalized(q);
        TRACEN("edge intersects facet on plane "<<plane(f)<<" on "<<q);
        call_back( e0, Object_handle(f), q);
      }
    }
  }
};

CGAL_END_NAMESPACE

#endif // SNC_POINT_LOCATOR_H

