// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// file          : include/CGAL/IO/Qt_widget_Nef_2.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================


#ifndef CGAL_QT_WIDGET_NEF_2_H
#define CGAL_QT_WIDGET_NEF_2_H

#include <CGAL/Nef_polyhedron_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/IO/Qt_widget.h>


typedef CGAL::Cartesian<double>::Point_2 Draw_point;
typedef CGAL::Cartesian<double>::Segment_2 Segment;
CGAL::Qt_widget	_W;
namespace CGAL{

template<typename T>
CGAL::Qt_widget& operator<<(CGAL::Qt_widget& ws, const Nef_polyhedron_2<T>& P)
{

}

void draw(Vertex_const_handle v) const
//{\Mop draws |v| according to the color and width specified by
//    |C.color(v)| and |C.width(v)|.}
{ 
	_W << CGAL::RED << point(v);
}

void draw(Halfedge_const_handle e) const
//{\Mop draws |e| according to the color and width specified by
//    |C.color(e)| and |C.width(e)|.}*/
{ 
  Segment s(point(source(e)),point(target(e));
  _W << CGAL::LineWidth(2);
  _W << CGAL::GREEN << s;
}

void get_point_list(std::list<Draw_point>& L, 
                    Halfedge_const_iterator e) const
{
  Halfedge_around_face_const_circulator fcirc(e), fend(fcirc);
  CGAL_For_all(fcirc,fend) {
    Point p = point(target(fcirc));
    L.push_back(Draw_point(CGAL::to_double(p.x()),
			   CGAL::to_double(p.y())));
  }   
}

void draw(Face_const_handle f) const
//{\Mop draws |f| with color |C.color(f)|.}
{ 
  //CGAL::Color cc = _CO.color(f,mark(f));
  //leda_color c (cc.r(),cc.g(),cc.b());
  std::list<Draw_point> outer_cycle;
  // First the outer face cycle:
  get_point_list(outer_cycle,halfedge(f));

  double x0 = _W.xmin();
  double y0 = _W.ymin();
  double x1 = _W.xmax();
  double y1 = _W.ymax();

  if ( _W.is_buffering() ) {
    x0 = _W.xreal(0);
    y0 = _W.yreal(_W.height());
    x1 = _W.xreal(_W.width());
    y1 = _W.yreal(0);
  }

  _W.reset_clip_mask();
  std::list<Draw_point> frame;
  frame.push_back(Draw_point(x0,y0));
  frame.push_back(Draw_point(x1,y0));
  frame.push_back(Draw_point(x1,y1));
  frame.push_back(Draw_point(x0,y1));
  frame.reverse();
  draw_face_cycle(frame,0); 
  // enforcing transparent mode outside f and inside frame
  draw_face_cycle(outer_cycle,1);
  // drawing non-transparent outer face cycle

  Hole_const_iterator hole_it;
  for (hole_it = holes_begin(f); hole_it != holes_end(f); ++hole_it) {
    std::list<Draw_point> hole;
    get_point_list(hole,hole_it);
    draw_face_cycle(hole,0);
    // enforcing transparent mode for holes
  }
  Isolated_vertex_const_iterator iv_it;
  for (iv_it = isolated_vertices_begin(f); 
       iv_it != isolated_vertices_end(f); ++iv_it) {
    draw(iv_it); 
  }
  _W.draw_box(x0,y0,x1,y1,c);
  _W.reset_clip_mask();
}


void draw_map() const
/*{\Mop draw the whole plane map.}*/
{
  Face_const_iterator 
    fit = faces_begin(), fend = faces_end();
  // we don't draw the first face outside the box:
  for ( ++fit; fit != fend; ++fit) 
    draw(fit);

  // draw segments underlying halfedges: 
  Halfedge_const_iterator hit, hend = halfedges_end();
  for (hit = halfedges_begin(); hit != hend; ++(++hit)) 
    draw(hit);

  // draw points underlying vertices:
  Vertex_const_iterator vit, vend = vertices_end();
  for (vit = vertices_begin(); vit != vend; ++vit) 
    draw(vit);
}

}//end namespace CGAL

#endif
