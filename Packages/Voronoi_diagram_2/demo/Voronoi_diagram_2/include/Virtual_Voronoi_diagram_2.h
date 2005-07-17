#ifndef CGAL_VIRTUAL_VORONOI_DIAGRAM_2_H
#define CGAL_VIRTUAL_VORONOI_DIAGRAM_2_H 1

#include <CGAL/basic.h>
#include <CGAL/iterator.h>
#include "../typedefs.h"

#ifdef CGAL_USE_QT
#include <CGAL/IO/Qt_widget.h>
#endif


CGAL_BEGIN_NAMESPACE


struct Virtual_Voronoi_diagram_2
{
  typedef CGAL::Object     Object;
  typedef ::Rep::Point_2   Point_2;
  typedef ::Rep::Circle_2  Circle_2;

  // insert a site
  virtual void insert(const Point_2&) = 0;
  virtual void insert(const Circle_2&) = 0;

#ifdef CGAL_USE_QT
  virtual void draw_feature(const Object&, Qt_widget&) const = 0;
  virtual void draw_diagram(Qt_widget&) const = 0;
  virtual void draw_sites(Qt_widget&) const = 0;
#endif  

  virtual Object locate(const Point_2&) const = 0;

  virtual Object ptr() = 0;

  virtual void clear() = 0;
};

//=========================================================================

template<class VD>
class Virtual_Voronoi_diagram_base_2
  : public VD, public Virtual_Voronoi_diagram_2
{
 protected:
  typedef Virtual_Voronoi_diagram_2    VBase;
  typedef VD                           Base;

  typedef typename VBase::Object       Object;
  typedef typename VBase::Point_2      Point_2;
  typedef typename VBase::Circle_2     Circle_2;

  Virtual_Voronoi_diagram_base_2() {}
  virtual ~Virtual_Voronoi_diagram_base_2() {}

  virtual void insert(const Point_2&) {}
  virtual void insert(const Circle_2&) {}

 public:
#ifdef CGAL_USE_QT
  virtual void draw_feature(const Object& o, Qt_widget& widget) const {
    typename Base::Locate_result lr;
    if ( !assign(lr, o) ) { return; }

    if ( lr.is_face() ) {
      typename Base::Face_handle f = lr;
      typename Base::Ccb_halfedge_circulator ccb_start = f->outer_ccb();
      typename Base::Ccb_halfedge_circulator ccb = ccb_start;
      do {
	typename Base::Voronoi_traits::Voronoi_edge_2 ve = ccb->curve();
	widget << ve;
	++ccb;
      } while ( ccb != ccb_start );
    }
  }

  virtual void draw_sites(Qt_widget& widget) const
  {
    for (typename Base::Generator_iterator git = this->generators_begin();
	 git != this->generators_end(); ++git) {
      widget << *git;
    }
  }

  virtual void draw_diagram(Qt_widget& widget) const
  {
    typename Base::Edge_iterator it;
    for (it = this->edges_begin(); it != this->edges_end(); ++it) {
      typename Base::Voronoi_traits::Voronoi_edge_2 ve2 = it->curve();
      widget << ve2;
    }
  }
#endif  

  virtual Object locate(const Point_2& q) const {
    typename Base::Voronoi_traits::Point_2 p(q.x(), q.y());
    typename Base::Locate_result lr = Base::locate(p);
    return CGAL::make_object(lr);
  }

  virtual Object ptr() = 0;

  virtual void clear() {
    Base::clear();
  }
};

//=========================================================================

class Concrete_Voronoi_diagram_2
  : public Virtual_Voronoi_diagram_base_2<VD2>
{
 protected:
  typedef Virtual_Voronoi_diagram_base_2<VD2> VBase;

  typedef VBase::Object   Object;
  typedef VBase::Base     Base;
  typedef VBase::Point_2  Point_2;

 public:
  Concrete_Voronoi_diagram_2() {}
  virtual ~Concrete_Voronoi_diagram_2() {}

  virtual void insert(const Point_2& p) {
    Base::Point_2 pp(p.x(), p.y());
    Base::insert(pp);
  }

  virtual Object ptr() { return CGAL::make_object(this); }
};

//=========================================================================

class Concrete_power_diagram_2
  : public Virtual_Voronoi_diagram_base_2<PD2>
{
 protected:
  typedef Virtual_Voronoi_diagram_base_2<PD2> VBase;

  typedef VBase::Object    Object;
  typedef VBase::Base      Base;
  typedef VBase::Point_2   Point_2;
  typedef VBase::Circle_2  Circle_2;

 public:
  Concrete_power_diagram_2() {}
  virtual ~Concrete_power_diagram_2() {}

  virtual void insert(const Point_2& p) {
    Base::Point_2 pp(p.x(), p.y());
    Base::Geom_traits::Weighted_point_2 wp(pp, 0);
    Base::insert(wp);
  }

  virtual void insert(const Circle_2& c) {
    Point_2 center = c.center();
    Base::Point_2 p(center.x(), center.y());
    Base::Geom_traits::Weighted_point_2 wp(p, c.squared_radius());
    Base::insert(wp);
  }

#ifdef CGAL_USE_QT
  virtual void draw_sites(Qt_widget& widget) const
  {
    VBase::draw_sites(widget);

    Base::Delaunay_graph::Finite_vertices_iterator vit;
    for (vit = this->dual().finite_vertices_begin();
	 vit != this->dual().finite_vertices_end(); ++vit) {
      Base::Geom_traits::Weighted_point_2 wp = vit->point();
      Point_2 center( CGAL::to_double(wp.point().x()),
		      CGAL::to_double(wp.point().y()) );
      Circle_2 c( center, CGAL::to_double(wp.weight()) );
      widget << c;
    }
  }
#endif

  virtual Object ptr() { return CGAL::make_object(this); }
};

//=========================================================================

class Concrete_Apollonius_diagram_2
  : public Virtual_Voronoi_diagram_base_2<AD2>
{
 protected:
  typedef Virtual_Voronoi_diagram_base_2<AD2> VBase;

  typedef VBase::Object    Object;
  typedef VBase::Base      Base;
  typedef VBase::Point_2   Point_2;
  typedef VBase::Circle_2  Circle_2;

 public:
  Concrete_Apollonius_diagram_2() {}
  virtual ~Concrete_Apollonius_diagram_2() {}

  virtual void insert(const Point_2& p) {
    Base::Geom_traits::Point_2 pp(p.x(), p.y());
    Base::Geom_traits::Site_2 s(p, 0);
    Base::insert(s);
  }

  virtual void insert(const Circle_2& c) {
    ::Rep::Point_2 center = c.center();
    Base::Geom_traits::Point_2 p(center.x(), center.y());
    Base::Geom_traits::Site_2::Weight w = CGAL::sqrt(c.squared_radius());
    Base::Geom_traits::Site_2 s(p, w);
    Base::insert(s);
  }

  virtual Object ptr() { return CGAL::make_object(this); }
};


CGAL_END_NAMESPACE


#endif // CGAL_VIRTUAL_VORONOI_DIAGRAM_2_H
