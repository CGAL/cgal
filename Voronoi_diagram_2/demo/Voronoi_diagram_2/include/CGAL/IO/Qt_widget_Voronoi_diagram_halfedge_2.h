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

#ifndef CGAL_QT_WIDGET_VORONOI_DIAGRAM_HALFEDGE_2_H
#define CGAL_QT_WIDGET_VORONOI_DIAGRAM_HALFEDGE_2_H 1


#include <CGAL/IO/Qt_widget.h>

namespace CGAL {

template<class VDA>
class Voronoi_diagram_halfedge_2
  : public VDA::Halfedge
{
 protected:
  typedef VDA                                        Voronoi_diagram;
  typedef typename Voronoi_diagram::Delaunay_graph   Delaunay_triangulation_2;
  typedef typename Voronoi_diagram::Halfedge         Base;
  typedef typename Base::Delaunay_edge               Delaunay_edge;

  typedef typename Voronoi_diagram::Adaptation_traits::Site_2  Site_2;

 public:
  Voronoi_diagram_halfedge_2() : Base() {}
  Voronoi_diagram_halfedge_2(const Base& e) : Base(e), is_conflict(false) {}
  Voronoi_diagram_halfedge_2(const Delaunay_edge& e, int inf, const Site_2& s)
    : Base(), is_conflict(true), e_(e), inf_(inf), s_(s) {}

  void draw(Qt_widget& qt_w) const
  {
    typedef typename Delaunay_triangulation_2::Geom_traits   Geom_traits;
    typedef typename Geom_traits::Point_2                    Point_2;
    typedef typename Geom_traits::Segment_2                  Segment_2;
    typedef typename Geom_traits::Line_2                     Line_2;

    if ( is_conflict ) {
      if ( inf_ == 0 ) {
	typename Geom_traits::Construct_circumcenter_2 circumcenter;
	Point_2 c1 = circumcenter(e_.first->vertex(0)->point(),
				  e_.first->vertex(1)->point(),
				  e_.first->vertex(2)->point());
	int ccw_i = (e_.second + 1) % 3;
	int cw_i  = (e_.second + 2) % 3;
	Point_2 c2 = circumcenter(e_.first->vertex(ccw_i)->point(),
				  e_.first->vertex(cw_i)->point(),
				  s_);
	qt_w << Segment_2(c1, c2);
      } else {
	typename Geom_traits::Construct_circumcenter_2 circumcenter;
	typename Geom_traits::Construct_bisector_2     c_bis;
	typename Geom_traits::Construct_ray_2          c_ray;
	int ccw_i = (e_.second + 1) % 3;
	int cw_i  = (e_.second + 2) % 3;
	Point_2 c = circumcenter(e_.first->vertex(ccw_i)->point(),
				 e_.first->vertex(cw_i)->point(),
				 s_);
	Line_2 l = c_bis(e_.first->vertex(ccw_i)->point(),
			 e_.first->vertex(cw_i)->point());
	qt_w << c_ray(c, l);
      }
      return;
    }

    if ( this->has_source() && this->has_target() ) {
      typename Geom_traits::Construct_circumcenter_2 circumcenter;
      Point_2 c1 = circumcenter(this->down()->point(),
				this->up()->point(),
				this->left()->point());
      Point_2 c2 = circumcenter(this->up()->point(),
				this->down()->point(),
				this->right()->point());
      qt_w << Segment_2(c1, c2);
    } else if ( this->has_source() && !this->has_target() ) {
      typename Geom_traits::Construct_circumcenter_2 circumcenter;
      typename Geom_traits::Construct_bisector_2     c_bis;
      typename Geom_traits::Construct_ray_2          c_ray;
      Point_2 c = circumcenter(this->down()->point(),
			       this->up()->point(),
			       this->left()->point());
      Line_2 l = c_bis(this->up()->point(), this->down()->point());
      qt_w << c_ray(c, l);
    } else if ( !this->has_source() && this->has_target() ) {
      typename Geom_traits::Construct_circumcenter_2 circumcenter;
      typename Geom_traits::Construct_bisector_2     c_bis;
      typename Geom_traits::Construct_ray_2          c_ray;
      Point_2 c = circumcenter(this->up()->point(),
			       this->down()->point(),
			       this->right()->point());
      Line_2 l = c_bis(this->down()->point(), this->up()->point());
      qt_w << c_ray(c, l);
    } else {
      CGAL_assertion( !this->has_source() && !this->has_target() );
      typename Geom_traits::Construct_bisector_2 c_bis;
      qt_w << c_bis(this->up()->point(),
		    this->down()->point());
    }
  }

private:
  bool is_conflict;
  Delaunay_edge e_;
  int inf_;
  Site_2 s_;
};

template<class VDA>
Qt_widget& operator<<(Qt_widget& qt_w,
		      const Voronoi_diagram_halfedge_2<VDA>& e)
{
  e.draw(qt_w);
  return qt_w;
}


} //namespace CGAL


#endif // CGAL_QT_WIDGET_VORONOI_DIAGRAM_HALFEDGE_2_H
