#ifndef CGAL_SNC_POINT_LOCATOR_H
#define CGAL_SNC_POINT_LOCATOR_H

#include <CGAL/basic.h>
#include <CGAL/intersection_3.h>
#include <CGAL/predicates_on_points_3.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_constructor.h>
#include <CGAL/Nef_3/SNC_FM_decorator.h>
#include <CGAL/Nef_3/SNC_SM_point_locator.h>

CGAL_BEGIN_NAMESPACE

template <typename SNC_>
class SNC_point_locator : public SNC_decorator<SNC_>
{ 
public:
  typedef SNC_                             SNC_structure;
  typedef SNC_point_locator<SNC_>          Self;
  typedef CGAL::SNC_decorator<SNC_>        Base;
  typedef CGAL::SNC_constructor<SNC_>      SNC_constructor;
  typedef CGAL::SNC_FM_decorator<SNC_>     FM_decorator;
  typedef CGAL::SNC_SM_decorator<SNC_>     SM_decorator;
  typedef CGAL::SNC_SM_point_locator<SNC_> SM_point_locator;

#define USING(t) typedef typename Base::t t
  USING(Vertex_iterator);   USING(Vertex_handle);
  USING(Halfedge_iterator); USING(Halfedge_handle);
  USING(Halffacet_iterator);    USING(Halffacet_handle);
  USING(Volume_iterator);   USING(Volume_handle);
  USING(SVertex_iterator);   USING(SVertex_handle);
  USING(SHalfedge_iterator); USING(SHalfedge_handle);
  USING(SHalfloop_iterator); USING(SHalfloop_handle);
  USING(SFace_iterator);     USING(SFace_handle);
  USING(Object_handle);
  USING(SObject_handle);
  USING(SFace_cycle_iterator);
  USING(Halffacet_cycle_iterator);
  USING(Shell_entry_iterator);
  USING(SHalfedge_around_facet_circulator);
  USING(Point_3);
  USING(Segment_3);
  USING(Plane_3);
  USING(Vector_3);
  USING(Line_3);
  USING(Sphere_point);
  USING(Sphere_segment);
  USING(Sphere_circle);
  USING(Mark);
#undef USING
  typedef void* GenPtr;

  SNC_point_locator(SNC_structure& W) : Base(W) {}

  typename SNC_::Object_handle locate(const Point_3& p) const;
  typename SNC_::Halffacet_handle get_facet_below(const Vertex_handle& vi) 
    const;
  
  void shorten(Segment_3& s, const Point_3& p) const
  { s = Segment_3(s.source(),p); }

  bool contains_internally(const Segment_3& s, const Point_3& p) const
  { CGAL_assertion(s.has_on(p));
    CGAL::Comparison_result r1 = compare_xyz(s.source(),p); 
    CGAL::Comparison_result r2 = compare_xyz(s.target(),p); 
    return (r1 == opposite(r2));
  }

  bool do_intersect(const Segment_3& s,
		    Halfedge_handle e,
		    Point_3& p) const
  { CGAL::Object o = CGAL::intersection(Line_3(s),Line_3(segment(e))); 
    if ( !CGAL::assign(p,o) ) return false;
    CGAL::Oriented_side os1 = h.oriented_side(s.source());
    CGAL::Oriented_side os2 = h.oriented_side(s.target());
    return (os1 == CGAL::opposite(os2));
  }

  bool do_intersect(const Segment_3& s,
		    Halffacet_handle f,
		    Point_3& p) const
  { Plane_3 h(plane(f));
    CGAL::Object o = CGAL::intersection(h, Line_3(s));
    if ( !CGAL::assign(p,o) ) return false;
    CGAL::Oriented_side os1 = h.oriented_side(s.source());
    CGAL::Oriented_side os2 = h.oriented_side(s.target());
    return (os1 == CGAL::opposite(os2));
  }


}; // SNC_point_locator<SNC_>



template <typename SNC_>
typename SNC_::Object_handle
SNC_point_locator<SNC_>::locate(const Point_3& p) const
{
  if ( number_of_vertices() == 0 ) return Object_handle();
  // now not empty:

  Vertex_iterator v;
  CGAL_forall_vertices(v,*sncp()) 
    if ( p == point(v) ) return Object_handle(v);

  Halfedge_iterator e;
  CGAL_forall_edges(e,*sncp()) 
    if ( segment(e).constains(p) ) return Object_handle(e);

  Halffacet_iterator f;
  CGAL_forall_halffacets(f,*sncp()) { // nur ufacets
    FM_decorator D(f);
    if ( D.facet_contains(p) ) return Object_handle(f);
  }
  // now not contained in any boundary object:

  Segment_3 s(p, point(vertices_begin())), so(s.opposite());
  Object_handle closest(vertices_begin());
  CGAL_forall_vertices(v,*sncp())
    if ( s.has_on(point(v)) ) 
    { shorten(s,point(v)); closest=Object_handle(v); }
  CGAL_forall_edges(e,*sncp()) {
    Point_3 q;
    if ( do_intersect(s,e,q) ) 
    { shorten(s,q); closest=Object_handle(e); }
  }
  CGAL_forall_halffacets(f,*sncp()) {
    Point_3 q;
    if ( do_intersect(s,f,q) )
    { shorten(s,q); closest=Object_handle(f); }
  }

  // now closest holds the closest boundary object that
  // contains the point s.target() :

  if ( CGAL::assign(v,closest) ) {
    SM_point_locator L(v);
    Sphere_point pl(so.direction());
    SObject_handle o = L.locate(pl);
    SFace_handle f;  
    if ( !CGAL::assign(f,o) ) CGAL_nef3_assertion_msg(0,"Must be sface.");
    /* v is the closest object and therefore the view towards p must
       be an sface */
    return Object_handle(volume(f));
  } else if ( CGAL::assign(e,closest) ) {
    SNC_constructor C(*sncp());
    Vertex_handle v_tmp = C.create_from_edge(e,s.target());
    SM_point_locator L(v_tmp);
    Sphere_point pl(so.direction());
    SObject_handle o = L.locate(pl);
    SFace_handle f;  
    if ( !CGAL::assign(f,o) ) CGAL_nef3_assertion_msg(0,"Must be sface.");
    /* e is the closest object and therefore the view towards p must
       be an sface */
    return Object_handle(volume(f));
  } else if ( CGAL::assign(f,closest) ) {
    // sideness correction necessary:
    if ( plane(f).has_on_negative_side(p) ) f = twin(f);
    CGAL_nef3_assertion( plane(f).has_on_positive_side(p) );
    return Object_handle(volume(f));
  } else CGAL_nef3_assertion_msg(0,"damn generic handles.");
  return Object_handle(); // never reached
}

CGAL_END_NAMESPACE
#endif //CGAL_SNC_POINT_LOCATOR_H

