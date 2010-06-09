// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
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
// $URL$
// $Id$
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_QT_WIDGET_APOLLONIUS_DIAGRAM_HALFEDGE_2_H
#define CGAL_QT_WIDGET_APOLLONIUS_DIAGRAM_HALFEDGE_2_H 1


#include <CGAL/IO/Qt_widget.h>
#include <CGAL/Apollonius_graph_2/Constructions_C2.h>
#include <CGAL/Hyperbola_segment_2.h>
#include <CGAL/Hyperbola_ray_2.h>
#include <CGAL/Hyperbola_2.h>

namespace CGAL {

template<class VDA>
class Apollonius_diagram_halfedge_2
  : public VDA::Halfedge
{
 protected:
  typedef VDA                                        Voronoi_diagram;
  typedef typename Voronoi_diagram::Delaunay_graph   Apollonius_graph_2;
  typedef typename Voronoi_diagram::Halfedge         Base;
  typedef typename Base::Delaunay_edge               Delaunay_edge;
  typedef typename Apollonius_graph_2::Geom_traits   Geom_traits;

  typedef typename Voronoi_diagram::Adaptation_traits::Site_2   Site_2;
  typedef typename Voronoi_diagram::Adaptation_traits::Point_2  Point_2;


  typedef CGAL_APOLLONIUS_GRAPH_2_NS::
  Construct_Apollonius_bisector_2<Geom_traits>
  Construct_Apollonius_bisector_2;

  typedef CGAL_APOLLONIUS_GRAPH_2_NS::
  Construct_Apollonius_bisector_ray_2<Geom_traits>
  Construct_Apollonius_bisector_ray_2;

  typedef CGAL_APOLLONIUS_GRAPH_2_NS::
  Construct_Apollonius_bisector_segment_2<Geom_traits>
  Construct_Apollonius_bisector_segment_2;


 public:
  Apollonius_diagram_halfedge_2() : Base() {}
  Apollonius_diagram_halfedge_2(const Base& e)
    : Base(e), is_conflict(false) {}
  Apollonius_diagram_halfedge_2(const Delaunay_edge& e, int inf,
				const Site_2& s)
    : Base(), is_conflict(true), e_(e), inf_(inf), s_(s) {}

  void draw(Qt_widget& qt_w) const
  {
    typedef typename Geom_traits::Assign_2               Assign_2;
    typedef typename Geom_traits::Segment_2              Segment_2;
    typedef typename Geom_traits::Ray_2                  Ray_2;
    typedef typename Geom_traits::Line_2                 Line_2;
    typedef Hyperbola_segment_2<Geom_traits>             Hyperbola_segment_2;
    typedef Hyperbola_ray_2<Geom_traits>                 Hyperbola_ray_2;
    typedef Hyperbola_2<Geom_traits>                     Hyperbola_2;

    Assign_2 assign = Geom_traits().assign_2_object();
    Hyperbola_segment_2 hs;
    Hyperbola_ray_2 hr;
    Hyperbola_2 h;
    Segment_2 s;
    Ray_2 r;
    Line_2 l;

    Object o;
    if ( is_conflict ) {
      int ccw_i = (e_.second + 1) % 3;
      int cw_i  = (e_.second + 2) % 3;
      typename Geom_traits::Construct_Apollonius_vertex_2 cvertex =
	Geom_traits().construct_Apollonius_vertex_2_object();
      if ( inf_ == 0 ) {
	Point_2 c1 = cvertex(e_.first->vertex(ccw_i)->site(),
			     e_.first->vertex(cw_i)->site(),
			     s_);

	Point_2 c2 = cvertex(e_.first->vertex(ccw_i)->site(),
			     e_.first->vertex(cw_i)->site(),
			     e_.first->vertex(e_.second)->site());

	Construct_Apollonius_bisector_segment_2 c_seg;
	o = c_seg(e_.first->vertex(ccw_i)->site(),
		  e_.first->vertex(cw_i)->site(),
		  c1, c2);
      } else if ( inf_ == 1 ) {
	Point_2 c = cvertex(e_.first->vertex(ccw_i)->site(),
			    e_.first->vertex(cw_i)->site(),
			    s_);

	Construct_Apollonius_bisector_ray_2 c_ray;
	o = c_ray(e_.first->vertex(ccw_i)->site(),
		  e_.first->vertex(cw_i)->site(),
		  c, POSITIVE);
      } else {
	CGAL_assertion( inf_ == 2 );
	Point_2 c1 = cvertex(e_.first->vertex(ccw_i)->site(),
			     e_.first->vertex(cw_i)->site(),
			     s_);

	Point_2 c2 = cvertex(e_.first->vertex(cw_i)->site(),
			     e_.first->vertex(ccw_i)->site(),
			     s_);

	Construct_Apollonius_bisector_segment_2 c_seg;
	o = c_seg(e_.first->vertex(ccw_i)->site(),
		  e_.first->vertex(cw_i)->site(),
		  c1, c2);
      }

      // fix this and use the output operators...
      if      ( assign(hs,o) )   hs.draw(qt_w);
      else if ( assign(hr,o) )   hr.draw(qt_w);
      else if ( assign(h, o) )   h.draw(qt_w);
      else if ( assign(s, o) )   qt_w << s;
      else if ( assign(r, o) )   qt_w << r;
      else if ( assign(l, o) )   qt_w << l;
      return;
    }

    if ( this->has_source() && this->has_target() ) {
      Construct_Apollonius_bisector_segment_2 c_seg;
      o = c_seg(this->down()->site(),
		this->up()->site(),
		this->left()->site(),
		this->right()->site());
    } else if ( this->has_source() && !this->has_target() ) {
      Construct_Apollonius_bisector_ray_2 c_ray;
      o = c_ray(this->down()->site(),
		this->up()->site(),
		this->left()->site());
    } else if ( !this->has_source() && this->has_target() ) {
      Construct_Apollonius_bisector_ray_2 c_ray;
      o = c_ray(this->up()->site(),
		this->down()->site(),
		this->right()->site());
    } else {
      CGAL_assertion( !this->has_source() && !this->has_target() );
      Construct_Apollonius_bisector_2 c_bis;
      o = c_bis(this->up()->site(),
		this->down()->site());
    }

    // fix this and use the output operators...
    if      ( assign(hs,o) )   hs.draw(qt_w);
    else if ( assign(hr,o) )   hr.draw(qt_w);
    else if ( assign(h, o) )   h.draw(qt_w);
    else if ( assign(s, o) )   qt_w << s;
    else if ( assign(r, o) )   qt_w << r;
    else if ( assign(l, o) )   qt_w << l;
  }

private:
  bool is_conflict;
  Delaunay_edge e_;
  int inf_;
  Site_2 s_;
};

template<class VDA>
Qt_widget& operator<<(Qt_widget& qt_w,
		      const Apollonius_diagram_halfedge_2<VDA>& e)
{
  e.draw(qt_w);
  return qt_w;
}


} //namespace CGAL


#endif // CGAL_QT_WIDGET_APOLLONIUS_DIAGRAM_HALFEDGE_2_H
