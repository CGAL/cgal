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
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/IO/Qt_widget.h>

#include <qnamespace.h>

typedef CGAL::Cartesian<double>::Point_2 Draw_point;
typedef CGAL::Cartesian<double>::Segment_2 Segment;


namespace CGAL{

template <typename T>
void draw(CGAL::Qt_widget& ws, 
	  Nef_polyhedron_2<T>::Vertex_const_iterator v);

template <typename T>
CGAL::Qt_widget& operator<<(CGAL::Qt_widget& ws, const Nef_polyhedron_2<T>& P)
{
    typedef Nef_polyhedron_2<T> Polyhedron;
    typedef typename T::RT RT;
    typedef typename T::Standard_RT Standard_RT;
    typedef typename Polyhedron::Topological_explorer 
                                    TExplorer;
    typedef typename TExplorer::Halfedge_around_face_const_circulator 
      Halfedge_around_face_const_circulator;
    typedef typename TExplorer::Hole_const_iterator
      Hole_const_iterator;
    typedef typename TExplorer::Isolated_vertex_const_iterator
      Isolated_vertex_const_iterator;

    typedef typename T::Standard_point_2 Standard_point_2;
    typedef typename T::Standard_segment_2 Standard_segment_2;
    typedef typename T::Standard_ray_2 Standard_ray_2;

    typedef Nef_polyhedron_2<T>::Vertex_const_handle
      Vertex_const_handle;
    typedef Nef_polyhedron_2<T>::Halfedge_const_handle
      Halfedge_const_handle;
    typedef Nef_polyhedron_2<T>::Face_const_handle
      Face_const_handle;

    typedef Nef_polyhedron_2<T>::Vertex_const_iterator
      Vertex_const_iterator;
    typedef Nef_polyhedron_2<T>::Halfedge_const_iterator
      Halfedge_const_iterator;
    typedef Nef_polyhedron_2<T>::Face_const_iterator
      Face_const_iterator;

    T traits;
    TExplorer D = P.explorer();
    const T& E = Polyhedron::EK;

    Standard_RT frame_radius = 100;
    E.determine_frame_radius(D.points_begin(), D.points_end(),
			     frame_radius);
    RT::set_R(frame_radius);
    
    //The faces
    Face_const_iterator 
      fit = D.faces_begin(), fend = D.faces_end();
    // we don't draw the first face outside the box:
    for ( ++fit; fit != fend; ++fit) {
      if(D.mark(fit))
	ws << CGAL::GREEN;
      else
	ws << CGAL::RED;
      
      std::list<Standard_point_2> l;
      Halfedge_around_face_const_circulator fcirc(D.halfedge(fit)), 
                                            fend(fcirc);
      CGAL_For_all(fcirc, fend){
	if(traits.is_standard(D.point(D.target(fcirc))))
	   l.push_back(traits.standard_point(D.point(D.target(fcirc))));
      }
      QPointArray array(l.size());int i=0;
      std::list<Standard_point_2>::const_iterator it = l.begin();
      while(it!=l.end()){
	array.setPoint(i++, ws.x_pixel(to_double((*it).x())), 
		       ws.y_pixel(to_double((*it).y())));
	it++;
      }
      ws.get_painter().drawPolygon(array);

      Hole_const_iterator hole_it;
      for (hole_it = D.holes_begin(fit); 
	   hole_it != D.holes_end(fit); ++hole_it) {
	std::list<Standard_point_2> hole;
	Halfedge_around_face_const_circulator fcirc(hole_it), 
                                              fend(fcirc);
	CGAL_For_all(fcirc, fend){
	  if(traits.is_standard(D.point(D.target(fcirc))))
	      hole.push_back(traits.standard_point(D.point(D.target(fcirc))));
	}//end CGAL_For_all
	QPointArray array(hole.size());int i=0;
	std::list<Standard_point_2>::const_iterator it = hole.begin();
	while(it!=hole.end()){
	  array.setPoint(i++, ws.x_pixel(to_double((*it).x())), 
		       ws.y_pixel(to_double((*it).y())));
	  it++;
	};

	if(D.mark(hole_it))
	  ws << CGAL::GREEN;
	else
	  ws << CGAL::RED;
	Qt::RasterOp old = ws.rasterOp();
	ws.setRasterOp(Qt::XorROP);
	ws.get_painter().drawPolygon(array);
	ws.setRasterOp(old);
      }//endfor
      
      ws << CGAL::RED << PointSize(5) << PointStyle(PLUS);      
      Isolated_vertex_const_iterator iv_it;
      for(iv_it = D.isolated_vertices_begin(fit);
	  iv_it != D.isolated_vertices_end(fit); ++iv_it)
	ws << traits.standard_point(D.point(iv_it));
    }//endfor Face_const_iterator
    
    // draw segments underlying halfedges: 
    Halfedge_const_iterator hit, hend = D.halfedges_end();
    for (hit = D.halfedges_begin(); hit != hend; ++(++hit)) {
      if(D.mark(hit))
	ws << CGAL::WHITE;
      else
	ws << CGAL::GRAY;
      if(traits.is_standard(D.point(D.source(hit))) 
	 && traits.is_standard(D.point(D.source(hit))))
      ws << Standard_segment_2(
		 traits.standard_point(D.point(D.source(hit))),
                 traits.standard_point(D.point(D.target(hit))));
    }
    
    // draw points underlying vertices:
    ws << CGAL::YELLOW << PointStyle(DISC);
    Vertex_const_iterator vit, vend = D.vertices_end();
    for (vit = D.vertices_begin(); vit != vend; ++vit){
      if(D.mark(vit))
	ws << CGAL::YELLOW;
      else
	ws << CGAL::GRAY;
      if(traits.is_standard(D.point(vit)))
            ws << traits.standard_point(D.point(vit));
    }
    
    return ws;
}

}//end namespace CGAL

#endif
