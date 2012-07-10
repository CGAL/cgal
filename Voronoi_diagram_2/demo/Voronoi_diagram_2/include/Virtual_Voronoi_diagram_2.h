// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_VIRTUAL_VORONOI_DIAGRAM_2_H
#define CGAL_VIRTUAL_VORONOI_DIAGRAM_2_H 1

#include <CGAL/basic.h>
#include <CGAL/iterator.h>
#include "../typedefs.h"

#include <CGAL/IO/Qt_widget.h>

#include <CGAL/IO/Qt_widget_Apollonius_diagram_halfedge_2.h>
#include <CGAL/IO/Qt_widget_Voronoi_diagram_halfedge_2.h>
#include <CGAL/IO/Qt_widget_power_diagram_halfedge_2.h>

namespace CGAL {


struct Virtual_Voronoi_diagram_2
{
  typedef CGAL::Object     Object;
  typedef ::Rep GlobalRep;
  typedef GlobalRep ::Point_2   Point_2;
  typedef GlobalRep::Circle_2  Circle_2;

  // virtual destructor
  virtual ~Virtual_Voronoi_diagram_2() {}

  // insert a site
  virtual void insert(const Point_2&) = 0;
  virtual void insert(const Circle_2&) = 0;

  // remove a site
  virtual void remove(const Object&) = 0;

  virtual void draw_feature(const Object&, Qt_widget&) const = 0;
  virtual void draw_diagram(Qt_widget&) const = 0;
  virtual void draw_sites(Qt_widget&) const = 0;
  virtual void draw_conflicts(const Point_2&, const Object&,
			      Qt_widget&) const = 0;
  virtual void draw_conflicts(const Circle_2&, const Object&,
			      Qt_widget&) const = 0;

  virtual Object locate(const Point_2&) const = 0;

  virtual Object get_conflicts(const Point_2&) const = 0;
  virtual Object get_conflicts(const Circle_2&) const = 0;

  virtual Object ptr() = 0;

  virtual bool is_valid() const = 0;

  virtual void clear() = 0;

};


//=========================================================================

template<class VD, class Halfedge_with_draw_t>
class Virtual_Voronoi_diagram_base_2
  : public VD, public Virtual_Voronoi_diagram_2
{
 protected:
  typedef Virtual_Voronoi_diagram_2    VBase;
  typedef VD                           Base;

  typedef typename VBase::Object       Object;
  typedef typename VBase::Point_2      Point_2;
  typedef typename VBase::Circle_2     Circle_2;

  typedef typename Base::Halfedge                 Halfedge;
  typedef typename Base::Halfedge_handle          Halfedge_handle;
  typedef typename Base::Face_handle              Face_handle;
  typedef typename Base::Ccb_halfedge_circulator  Ccb_halfedge_circulator;
  typedef typename Base::Edge_iterator            Edge_iterator;
  typedef typename Base::Site_iterator            Site_iterator;

  typedef typename Base::Locate_result            Locate_result;

  typedef Halfedge_with_draw_t                    Halfedge_with_draw;

  typedef typename Base::Delaunay_graph           Delaunay_graph;
  typedef typename Delaunay_graph::Edge           Delaunay_edge;
  typedef typename Delaunay_graph::Vertex_handle  Delaunay_vertex_handle;
  typedef typename Delaunay_graph::Face_handle    Delaunay_face_handle;

  typedef typename Base::Adaptation_traits::Site_2   Site_2;

  typedef Triangulation_cw_ccw_2                  CW_CCW_2;

  Virtual_Voronoi_diagram_base_2() {}
  virtual ~Virtual_Voronoi_diagram_base_2() {}

  virtual void insert(const Point_2&) {}
  virtual void insert(const Circle_2&) {}

  virtual void remove(const Object& /* o */) {
    // remove is not ready yet
#if 0
    Locate_result lr;
    if ( CGAL::assign(lr, o) && lr.is_face() ) {
      Face_handle f = lr;
      Base::remove(f);
    }
#endif
  }

  virtual
  Object conflicts(const Site_2& s) const
  {
    if ( Base::dual().dimension() < 2 ) {
      return CGAL::make_object( (int)0 );
    }

    typedef std::vector<Delaunay_edge>          Edge_vector;
    typedef std::vector<Delaunay_face_handle>   Face_vector;

    typedef std::back_insert_iterator<Face_vector>   Face_output_iterator;
    typedef std::back_insert_iterator<Edge_vector>   Edge_output_iterator;

    Edge_vector evec;
    Face_vector fvec;

    Face_output_iterator fit(fvec);
    Edge_output_iterator eit(evec);

    Base::dual().get_conflicts_and_boundary(s, fit, eit);

    return CGAL::make_object( std::make_pair(fvec, evec) );
  }

  Delaunay_edge opposite(const Delaunay_edge& e) const {
    int j = Base::dual().tds().mirror_index(e.first, e.second);
    Delaunay_face_handle n = e.first->neighbor(e.second);
    return Delaunay_edge(n, j);
  }

  template<class Query, class Iterator>
  bool find(const Query& q, Iterator first, Iterator beyond) const {
    for (Iterator it = first; it != beyond; ++it) {
      if ( q == *it ) { return true; }
    }
    return false;
  }

  virtual void draw_conflicts(const Site_2& s, const Object& o,
			      Qt_widget& widget) const {
    typedef std::vector<Delaunay_edge>           Edge_vector;
    typedef std::vector<Delaunay_face_handle>    Face_vector;
    typedef std::pair<Face_vector,Edge_vector>   result_type;

    result_type res;
    if ( !CGAL::assign(res, o) ) { return; }

    Face_vector fvec = res.first;
    Edge_vector evec = res.second;

    widget << CGAL::YELLOW;
    unsigned int linewidth = widget.lineWidth();
    widget << CGAL::LineWidth(4);

    bool do_regular_draw = true;
    if ( evec.size() == 2 ) {
      Delaunay_edge e1 = evec[0];
      Delaunay_edge e2 = evec[1];
      if ( e1 == opposite(e2) ) {
	do_regular_draw = false;
	if ( !Base::dual().is_infinite(e1) ) {
	  Halfedge_with_draw ee(e1, 2, s);
	  widget << ee;
	}
      }
    }

    if ( do_regular_draw ) {
      for (unsigned int i = 0; i < evec.size(); i++) {
	if ( Base::dual().is_infinite(evec[i]) ) { continue; }
	Delaunay_edge opp = opposite(evec[i]);
	Halfedge_with_draw ee(opp, Base::dual().is_infinite(opp.first), s);
	widget << ee;
      }
    }

    typename Base::Adaptation_policy::Edge_rejector e_rejector =
      Base::adaptation_policy().edge_rejector_object();
    for (unsigned int i = 0; i < fvec.size(); i++) {
      for (int j = 0; j < 3; j++) {
	Delaunay_edge e(fvec[i], j);
	Delaunay_edge opp = opposite(e);

	if ( Base::dual().is_infinite(e) ) { continue; }

	if ( !find(e, evec.begin(), evec.end()) &&
	     !find(opp, evec.begin(), evec.end()) ) {
	  if ( !e_rejector(Base::dual(),e) ) {
	    Halfedge_with_draw ee(*Base::dual(e));
	    widget << ee;
	  }
	}
      }
    }

    widget << CGAL::LineWidth(linewidth);
  }

 public:
  void draw_edge(const Halfedge& e, Qt_widget& widget) const {
    Halfedge_with_draw ee(e);
    widget << ee;
  }

  virtual void draw_feature(const Object& o, Qt_widget& widget) const {
    Locate_result lr;
    //    if ( !CGAL::assign(lr, o) ) { return; }
#if 1
    const Locate_result* lrp0 = CGAL::object_cast<Locate_result>(&o);
    Locate_result* lrp = const_cast<Locate_result*>(lrp0);
    if ( lrp == NULL ) { return; }
#else
    try {
      lr = CGAL::object_cast<Locate_result>(o);
    } catch ( CGAL::Bad_object_cast ) {
      return;
    }
#endif

    if ( Face_handle* f = boost::get<Face_handle>(lrp) ) {
      Ccb_halfedge_circulator ccb_start = (*f)->outer_ccb();
      Ccb_halfedge_circulator ccb = ccb_start;
      do {
	draw_edge(*ccb, widget);
	++ccb;
      } while ( ccb != ccb_start );
    }
  }

  virtual void draw_sites(Qt_widget& widget) const
  {
    for (Site_iterator sit = this->sites_begin();
	 sit != this->sites_end(); ++sit) {
      widget << *sit;
    }
  }

  virtual void draw_diagram(Qt_widget& widget) const
  {
    Edge_iterator it;
    for (it = this->edges_begin(); it != this->edges_end(); ++it) {
      draw_edge(*it, widget);
    }
  }

  virtual void draw_conflicts(const Point_2& ,	const Object& ,
			      Qt_widget& ) const {}

  virtual void draw_conflicts(const Circle_2& , const Object& ,
			      Qt_widget& ) const {}

  virtual Object locate(const Point_2& q) const {
    if ( Base::number_of_faces() == 0 ) {
      return CGAL::make_object(int(0));
    }
    typename Base::Adaptation_traits::Point_2 p(q.x(), q.y());
    Locate_result lr = Base::locate(p);
    return CGAL::make_object(lr);
  }

  virtual Object get_conflicts(const Point_2& ) const {
    return CGAL::make_object((int)0);
  }

  virtual Object get_conflicts(const Circle_2& ) const {
    return CGAL::make_object((int)0);
  }

  virtual Object ptr() = 0;

  virtual bool is_valid() const {
    return Base::is_valid();
  }

  virtual void clear() {
    Base::clear();
  }
};

//=========================================================================

class Concrete_Voronoi_diagram_2
  : public Virtual_Voronoi_diagram_base_2
  <VD2,Voronoi_diagram_halfedge_2<VD2> >
{
 protected:
  typedef Voronoi_diagram_halfedge_2<VD2>                   VD2_Halfedge;
  typedef Virtual_Voronoi_diagram_base_2<VD2,VD2_Halfedge>  VBase;

  typedef VBase::Object   Object;
  typedef VBase::Base     Base;
  typedef VBase::Point_2  Point_2;

  Base::Adaptation_traits::Site_2 to_site(const Point_2& p) const {
    return Base::Adaptation_traits::Site_2(p.x(), p.y());
  }

 public:
  Concrete_Voronoi_diagram_2() {}
  virtual ~Concrete_Voronoi_diagram_2() {}

  virtual void insert(const Point_2& p) {
    Base::insert( to_site(p) );
  }

  virtual Object get_conflicts(const Point_2& q) const {
    Base::Adaptation_traits::Point_2 p = to_site(q);
    return conflicts( to_site(q) );
  }

  virtual Object get_conflicts(const Circle_2& ) const {
    return CGAL::make_object( (int)0 );
  }

  virtual void draw_conflicts(const Point_2& p, const Object& o,
  			      Qt_widget& widget) const
  {
    VBase::draw_conflicts( to_site(p), o, widget);
  }

  virtual void draw_conflicts(const Circle_2& , const Object& ,
			      Qt_widget& ) const {
    return;
  }

  virtual Object ptr() { return CGAL::make_object(this); }
};

//=========================================================================

class Concrete_power_diagram_2
  : public Virtual_Voronoi_diagram_base_2
  <PD2,Power_diagram_halfedge_2<PD2> >
{
 protected:
  typedef Power_diagram_halfedge_2<PD2>                     PD2_Halfedge;
  typedef Virtual_Voronoi_diagram_base_2<PD2,PD2_Halfedge>  VBase;

  typedef VBase::Object    Object;
  typedef VBase::Base      Base;
  typedef VBase::Point_2   Point_2;
  typedef VBase::Circle_2  Circle_2;

  typedef Base::Delaunay_graph::Geom_traits  Geom_traits;

  Geom_traits::Weighted_point_2 to_site(const Point_2& p) const {
    Base::Point_2 pp(p.x(), p.y());
    return Geom_traits::Weighted_point_2(pp, 0);
  }

  Geom_traits::Weighted_point_2 to_site(const Circle_2& c) const {
    Point_2 center = c.center();
    Base::Point_2 p(center.x(), center.y());
    return Geom_traits::Weighted_point_2(p, c.squared_radius());
  }

 public:
  Concrete_power_diagram_2() {}
  virtual ~Concrete_power_diagram_2() {}

  virtual void insert(const Point_2& p) {
    Base::insert( to_site(p) );
  }

  virtual void insert(const Circle_2& c) {
    Base::insert( to_site(c) );
  }

  virtual Object get_conflicts(const Point_2& q) const {
    return conflicts( to_site(q) );
  }

  virtual Object get_conflicts(const Circle_2& c) const {
    return conflicts( to_site(c) );
  }

  virtual void draw_sites(Qt_widget& widget) const
  {
    VBase::draw_sites(widget);

    Base::Delaunay_graph::Finite_vertices_iterator vit;
    for (vit = this->dual().finite_vertices_begin();
	 vit != this->dual().finite_vertices_end(); ++vit) {
      Geom_traits::Weighted_point_2 wp = vit->point();
      Point_2 center( CGAL::to_double(wp.point().x()),
		      CGAL::to_double(wp.point().y()) );
      Circle_2 c( center, CGAL::to_double(wp.weight()) );
      widget << c;
    }
  }

  virtual void draw_conflicts(const Point_2& p, const Object& o,
  			      Qt_widget& widget) const
  {
    VBase::draw_conflicts( to_site(p), o, widget);
  }

  virtual void draw_conflicts(const Circle_2& c, const Object& o,
  			      Qt_widget& widget) const
  {
    VBase::draw_conflicts( to_site(c), o, widget);
  }

  virtual Object ptr() { return CGAL::make_object(this); }
};

//=========================================================================

class Concrete_Apollonius_diagram_2
  : public Virtual_Voronoi_diagram_base_2
  <AD2,Apollonius_diagram_halfedge_2<AD2> >
{
 protected:
  typedef Apollonius_diagram_halfedge_2<AD2>                AD2_Halfedge;
  typedef Virtual_Voronoi_diagram_base_2<AD2,AD2_Halfedge>  VBase;

  typedef VBase::Object    Object;
  typedef VBase::Base      Base;
  typedef VBase::Point_2   Point_2;
  typedef VBase::Circle_2  Circle_2;

  typedef Base::Delaunay_graph::Geom_traits  Geom_traits;

  Geom_traits::Site_2 to_site(const Point_2& p) const {
    Geom_traits::Point_2 pp(p.x(), p.y());
    return Geom_traits::Site_2(p, 0);
  }

  Geom_traits::Site_2 to_site(const Circle_2& c) const {
    typedef ::Rep GlobalRep;
    GlobalRep::Point_2 center = c.center();
    Geom_traits::Point_2 p(center.x(), center.y());
    Geom_traits::Site_2::Weight w = CGAL::sqrt(c.squared_radius());
    return Geom_traits::Site_2(p, w);
  }

 public:
  Concrete_Apollonius_diagram_2() {}
  virtual ~Concrete_Apollonius_diagram_2() {}

  virtual void insert(const Point_2& p) {
    Base::insert( to_site(p) );
  }

  virtual void insert(const Circle_2& c) {
    Base::insert( to_site(c) );
  }

  virtual Object get_conflicts(const Point_2& p) const {
    return conflicts( to_site(p) );
  }

  virtual Object get_conflicts(const Circle_2& c) const {
    return conflicts( to_site(c) );
  }

  virtual void draw_sites(Qt_widget& widget) const
  {
    VBase::draw_sites(widget);

    Base::Delaunay_graph::Hidden_sites_iterator hit;
    for (hit = this->dual().hidden_sites_begin();
	 hit != this->dual().hidden_sites_end(); ++hit) {
      widget << *hit;
    }
  }

  virtual void draw_conflicts(const Point_2& p, const Object& o,
  			      Qt_widget& widget) const
  {
    VBase::draw_conflicts( to_site(p), o, widget);
  }

  virtual void draw_conflicts(const Circle_2& c, const Object& o,
  			      Qt_widget& widget) const
  {
    VBase::draw_conflicts( to_site(c), o, widget);
  }

  virtual Object ptr() { return CGAL::make_object(this); }
};


} //namespace CGAL


#endif // CGAL_VIRTUAL_VORONOI_DIAGRAM_2_H
