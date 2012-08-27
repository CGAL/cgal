// Copyright (c) 2002-2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Radu Ursu

#ifndef CGAL_QT_WIDGET_NEF_2_H
#define CGAL_QT_WIDGET_NEF_2_H

#include <CGAL/Nef_polyhedron_2.h>
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/Quotient.h>

namespace CGAL{

template <class NT>
CGAL::Quotient<NT>
d_to_q(double x)
{ 
    NT num = 0; 
    NT den = 1;

    if (x != 0.0)
    { int neg = (x < 0);
      if (neg) x = -x;

      const unsigned shift = 15;   // a safe shift per step
      const unsigned int shift_pow = 32768; // = 2^shift
      const double width = 32768;  // = 2^shift
      const int maxiter = 20;      // ought not be necessary, but just in case,
                                   // max 300 bits of precision
      int expt;
      double mantissa = std::frexp(x, &expt);
      long exponent = expt;
      double intpart;
      int k = 0;
      
      while (mantissa != 0.0 && k++ < maxiter)

      { mantissa *= width; // shift double mantissa
        mantissa = std::modf(mantissa, &intpart);
        num *= (long)shift_pow;
        num += (long)intpart;
        exponent -= shift;
      }
      int expsign = (exponent>0 ? +1 : (exponent<0 ? -1 : 0));
      exponent *= expsign;
      NT twopot(2);
      NT exppot(1);
      while (exponent!=0) {
        if (exponent & 1)
          exppot *= twopot;
        exponent >>= 1;
        twopot *= twopot;
      }

      if (expsign > 0)
        num *= exppot;
      else if (expsign < 0)
        den *= exppot;
      if (neg)
        num = -num;
    }
    CGAL::Quotient<NT> q(num,den);
    q.normalize();
    return q;
}


template <typename T>
CGAL::Qt_widget& operator<<(CGAL::Qt_widget& ws, const Nef_polyhedron_2<T>& P)
{
    typedef Nef_polyhedron_2<T> Polyhedron;
    typedef typename T::Standard_RT Standard_RT;
    typedef typename T::Standard_segment_2
      Standard_segment_2;
    typedef typename T::Standard_point_2
      Standard_point_2;
    typedef typename Polyhedron::Explorer TExplorer;
    typedef typename TExplorer::Halfedge_around_face_const_circulator 
      Halfedge_around_face_const_circulator;

    typedef typename TExplorer::Point Point;

    typedef typename TExplorer::Vertex_const_iterator
      Vertex_const_iterator;
    typedef typename TExplorer::Halfedge_const_iterator
      Halfedge_const_iterator;
    typedef typename TExplorer::Face_const_iterator
      Face_const_iterator;

    //get the background color, fill color, and the object color
    QColor bgcolor = ws.backgroundColor();
	  QColor fillcolor = ws.fillColor();
	  QColor color = ws.color();

    //QPixmap
    //QPainter painter
    QPixmap &widget_pixmap = ws.get_pixmap();
    QPixmap copy_of_pixmap = (QPixmap)widget_pixmap;
    widget_pixmap.fill(bgcolor);


    //Get the screen rectangle to intersect with the current Nef
    CGAL::Quotient<Standard_RT> wsxq = d_to_q<Standard_RT>(ws.x_min()-1);
    CGAL::Quotient<Standard_RT> wsyq = d_to_q<Standard_RT>(ws.y_min()-1);
    Standard_RT wsx = wsxq.numerator() * wsyq.denominator(); 
    Standard_RT wsy = wsyq.numerator() * wsxq.denominator(); 
    Standard_RT wsh  = wsxq.denominator() * wsyq.denominator(); 
    Standard_point_2 p1(wsx, wsy, wsh);
    
    wsxq = d_to_q<Standard_RT>(ws.x_min()-1);
    wsyq = d_to_q<Standard_RT>(ws.y_max()+1);
    wsx = wsxq.numerator() * wsyq.denominator(); 
    wsy = wsyq.numerator() * wsxq.denominator(); 
    wsh  = wsxq.denominator() * wsyq.denominator(); 
    Standard_point_2 p2(wsx, wsy, wsh);

    wsxq = d_to_q<Standard_RT>(ws.x_max()+1);
    wsyq = d_to_q<Standard_RT>(ws.y_max()+1);
    wsx = wsxq.numerator() * wsyq.denominator(); 
    wsy = wsyq.numerator() * wsxq.denominator(); 
    wsh  = wsxq.denominator() * wsyq.denominator(); 
    Standard_point_2 p3(wsx, wsy, wsh);

    wsxq = d_to_q<Standard_RT>(ws.x_max()+1);
    wsyq = d_to_q<Standard_RT>(ws.y_min()-1);
    wsx = wsxq.numerator() * wsyq.denominator(); 
    wsy = wsyq.numerator() * wsxq.denominator(); 
    wsh  = wsxq.denominator() * wsyq.denominator(); 
    Standard_point_2 p4(wsx, wsy, wsh);

    Standard_point_2 rect1[4] = {p4, p3, p2, p1};
    Nef_polyhedron_2<T> N1(rect1, rect1+4);
    Nef_polyhedron_2<T> N2 = P.intersection(N1);
    TExplorer D = N2.explorer();

	  //TExplorer D = P.explorer();
    
    //The faces
    Face_const_iterator 
      fit = D.faces_begin(), fend = D.faces_end();
    // we don't draw the first face outside the box:
    for ( ++fit; fit != fend; ++fit) {
      Qt::RasterOp old_raster = ws.rasterOp();
      ws.setRasterOp(Qt::CopyROP);
      //save the initial raster mode
      if(D.mark(fit))
      	ws.setFillColor(fillcolor);      
      else
        ws.setFillColor(bgcolor);        

      std::list<Point> l;
      Halfedge_around_face_const_circulator fcirc(D.halfedge(fit)), 
                                            fend(fcirc);
      CGAL_For_all(fcirc, fend){
        if(D.is_standard(D.target(fcirc)))
        l.push_back(D.point(D.target(fcirc)));
      }
      QPointArray array(l.size());int i=0;
      typename std::list<Point>::const_iterator it = l.begin();
      while(it!=l.end()){
        array.setPoint(i++, ws.x_pixel(to_double((*it).x())),
		      ws.y_pixel(to_double((*it).y())));
      it++;
      }
      ws.get_painter().drawPolygon(array);
      ws.setRasterOp(old_raster);
/*
      typedef typename TExplorer::Isolated_vertex_const_iterator
      Isolated_vertex_const_iterator;
      Isolated_vertex_const_iterator iv_it;
      for (iv_it = D.isolated_vertices_begin(fit); 
        iv_it != D.isolated_vertices_end(fit); ++iv_it) {
        if(D.mark(iv_it))
          ws.setColor(color);
        else
       	  ws.setColor(bgcolor);
        if(D.is_standard(iv_it))
          ws << D.point(iv_it);
      }
*/
    }//endfor Face_const_iterator
    
    //save the initial raster mode
    Qt::RasterOp old_raster = ws.rasterOp();
    ws.setRasterOp(Qt::CopyROP);


    // draw segments underlying halfedges: 
    Halfedge_const_iterator hit, hend = D.halfedges_end();
    for (hit = D.halfedges_begin(); hit != hend; ++(++hit)) {
      if(D.mark(hit))
        ws.setColor(color);
      else
        ws.setColor(bgcolor);
      if(D.is_standard(D.source(hit)) 
        && D.is_standard(D.target(hit)))
          ws << Standard_segment_2(D.point(D.source(hit)),
                                   D.point(D.target(hit)));
    }
    
    // draw points underlying vertices:
    Vertex_const_iterator vit, vend = D.vertices_end();
    for (vit = D.vertices_begin(); vit != vend; ++vit){
      if(D.mark(vit))
        ws.setColor(color);
      else
       	ws.setColor(bgcolor);
      if(D.is_standard(vit))
        ws << D.point(vit);
    }
    
    ws.setRasterOp(old_raster);
    bitBlt(&widget_pixmap, 0, 0, &copy_of_pixmap, 
      0, 0, ws.width(), ws.height(), Qt::XorROP, true);
    return ws;
}

}//end namespace CGAL

#endif
