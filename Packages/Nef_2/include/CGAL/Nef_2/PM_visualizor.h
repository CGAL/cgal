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
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Nef_2/PM_visualizor.h
// package       : Nef_2 
// chapter       : Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>
//
// implementation: LEDA drawer
// ============================================================================

#ifndef PM_VISUALIZOR_H
#define PM_VISUALIZOR_H

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/IO/Window_stream.h>

#define USING(t) typedef typename PMCDEC::t t
#define LGREY CGAL::Color(190,190,190)
#define DGREY CGAL::Color(130,130,130)
CGAL_BEGIN_NAMESPACE

template <typename PMCDEC>
class PM_BooleColor 
{
  USING(Vertex_const_handle);   
  USING(Halfedge_const_handle); 
  USING(Face_const_handle);
  USING(Mark);
public:
  Color color(Vertex_const_handle, const Mark& m) const
  { return ( m ? CGAL::BLACK : LGREY ); }
  int width(Vertex_const_handle, const Mark& m) const
  { return 3; }
  Color color(Halfedge_const_handle, const Mark& m) const
  { return ( m ? CGAL::BLACK : LGREY ); }
  int width(Halfedge_const_handle, const Mark& m) const
  { return 2; }
  Color color(Face_const_handle, const Mark& m) const
  { return ( m ? DGREY : CGAL::WHITE ); }
};


template <typename PMCDEC>
class PM_DefColor 
{
  CGAL::Color _cs, _cf;
  int _wv, _we;
public:
  PM_DefColor() : 
    _cs(CGAL::BLACK), _cf(CGAL::WHITE), _wv(3), _we(2) {}
  PM_DefColor(CGAL::Color cs, CGAL::Color cf, int wv, int we) : 
    _cs(cs), _cf(cf), _wv(wv), _we(we) {}

  USING(Vertex_const_handle);   
  USING(Halfedge_const_handle); 
  USING(Face_const_handle);
  USING(Mark);
  Color color(Vertex_const_handle, const Mark&) const
  { return _cs; }
  int width(Vertex_const_handle, const Mark&) const
  { return _wv; }
  Color color(Halfedge_const_handle, const Mark&) const
  { return _cs; }
  int width(Halfedge_const_handle, const Mark&) const
  { return _we; }
  Color color(Face_const_handle, const Mark&) const
  { return _cf; }
};


/*{\Moptions outfile=PM_visualizor.man }*/
/*{\Manpage {PM_visualizor}{PMCDEC,GEOM,COLORDA}{Drawing plane maps}{V}}*/

template <typename PMCDEC, typename GEOM, 
          typename COLORDA = PM_DefColor<PMCDEC> >
class PM_visualizor : public PMCDEC
{
/*{\Mdefinition An instance |\Mvar| of the data type |\Mname| is a
decorator to draw the structure of a plane map into a CGAL window
stream. It is generic with respect to two template concepts.  |PMCDEC|
has to be a decorator model of our |PM_const_decorator|
concept. |GEOM| has to be a model of our geometry kernel concept.
The data accessor |COLORDA| has to have two members determining
the visualization parameters of the objects of |P|:\\
|CGAL::Color color(Vertex/Halfedge/Face_const_handle h) const|\\
|int width(Vertex/Halfedge_const_handle h) const|.
}*/

/*{\Mgeneralization PMCDEC}*/
/*{\Mtypes 3}*/
public:
  typedef PM_visualizor<PMCDEC,GEOM,COLORDA> Self;
  typedef PMCDEC Base;

  USING(Plane_map);
  USING(Vertex_const_handle);   
  USING(Halfedge_const_handle); 
  USING(Face_const_handle);     
  USING(Vertex_const_iterator);
  USING(Halfedge_const_iterator);
  USING(Face_const_iterator);
  USING(Halfedge_around_face_const_circulator);
  USING(Halfedge_around_vertex_const_circulator);
  USING(Hole_const_iterator);
  USING(Isolated_vertex_const_iterator);
  USING(Point);
  USING(Mark);
  typedef typename GEOM::Segment_2 Segment;
  typedef CGAL::Cartesian<double>::Point_2 Draw_point;

  typedef PMCDEC PM_const_decorator;
  /*{\Mtypemember The plane map decorator.}*/

  typedef GEOM Geometry;
  /*{\Mtypemember The used geometry.}*/

  typedef COLORDA Color_objects;
  /*{\Mtypemember The color data accessor.}*/

  CGAL::Window_stream& _W;
  const Geometry& _K;
  const Color_objects& _CO;

/*{\Mcreation 4}*/
PM_visualizor(CGAL::Window_stream& W, 
  const PM_const_decorator& D,
  const Geometry& K = Geometry(),
  const Color_objects& C = Color_objects() ) 
/*{\Mcreate creates an instance |\Mvar| of type |\Mname|
    to visualize the vertices, edges, and faces of |D| in window |W|.
    The coloring of the objects is determined by data accessor |C|.}*/
  : Base(D), _W(W), _K(K), _CO(C) 
{ _W.set_node_width(3); 
  _W.set_line_width(2);
}


/*{\Moperations 2 1}*/

void draw(Vertex_const_handle v) const
/*{\Mop draws |v| according to the color and width specified by
    |C.color(v)| and |C.width(v)|.}*/
{ int ow = _W.set_node_width(_CO.width(v,mark(v)));
  _W << _CO.color(v,mark(v)) << point(v);
  _W.set_node_width(ow);
}

void draw(Halfedge_const_handle e) const
/*{\Mop draws |e| according to the color and width specified by
    |C.color(e)| and |C.width(e)|.}*/
{ int ow = _W.set_line_width(_CO.width(e,mark(e)));
  Segment s = _K.construct_segment(point(source(e)),point(target(e)));
  _W << _CO.color(e,mark(e)) << s;
  _W.set_line_width(ow);
}

void draw_face_cycle(const std::list<Draw_point>& fc, int c) const
{ int n = fc.size();
  double* xc = new double[n];
  double* yc = new double[n];
  int i = 0;
  std::list<Draw_point>::const_iterator it;
  for (it = fc.begin(); it != fc.end(); ++i,++it) 
  { xc[i] = (*it).x(); yc[i] = (*it).y(); }
  _W.clip_mask_polygon(n,xc,yc,c);
  delete[] xc;
  delete[] yc;
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
/*{\Mop draws |f| with color |C.color(f)|.}*/
{ 
  CGAL::Color cc = _CO.color(f,mark(f));
  leda_color c (cc.r(),cc.g(),cc.b());
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

void init_window() const
{ _W.set_show_coordinates(true); _W.init(-110,110,-110); _W.display(); 
  _W.set_node_width(3); }

void draw_skeleton(const CGAL::Color& c=CGAL::BLACK) const
{ 
  int old = _W.set_line_width(1);
  _W << c;
  Halfedge_const_iterator hit, hend = halfedges_end();
  for (hit = halfedges_begin(); hit != hend; ++(++hit))
    _W << _K.construct_segment(point(source(hit)),point(target(hit)));

  // draw points underlying vertices:
  Vertex_const_iterator vit, vend = vertices_end();
  for (vit = vertices_begin(); vit != vend; ++vit) 
    _W << point(vit);
  _W.set_line_width(old);
}

void draw_ending_bundle(Vertex_const_handle v, 
        const CGAL::Color& c=CGAL::BLACK)
{ if (is_isolated(v)) return;
  _W << c;
  Halfedge_around_vertex_const_circulator hc(first_out_edge(v)), hend(hc);
  CGAL_For_all(hc,hend) {
    Point p1=point(source(hc)), p2=point(target(hc));
    if ( _K.compare_xy(p1,p2)>0 )
      _W << _K.construct_segment(p1,p2);
  }
}

}; // end of PM_visualizor 




CGAL_END_NAMESPACE
#undef USING
#undef LGREY
#undef DGREY
#endif // PM_VISUALIZOR_H

