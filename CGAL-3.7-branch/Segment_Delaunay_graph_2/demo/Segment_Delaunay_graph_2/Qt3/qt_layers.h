// Copyright (c) 2003,2004,2005  INRIA Sophia-Antipolis (France) and
// Notre Dame University (U.S.A.).  All rights reserved.
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
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>

#ifndef QT_LAYERS_H
#define QT_LAYERS_H

#include <CGAL/IO/Qt_widget_layer.h>

template< class T >
class Voronoi_diagram_layer : public CGAL::Qt_widget_layer {
private:
  T& sdg;

public:
  Voronoi_diagram_layer(T& sdg) : sdg(sdg) {}

  void draw() {
#if 1
    draw_diagram(*widget, sdg);
#else

    *widget << CGAL::BLUE;
#if !defined (__POWERPC__)
    *widget << CGAL::PointSize(3);
    *widget << CGAL::LineWidth(3);
#endif
    sdg.draw_dual(*widget);

#endif
  }

};

template< class T >
class Skeleton_layer : public CGAL::Qt_widget_layer {
private:
  T& sdg;

public:
  Skeleton_layer(T& sdg) : sdg(sdg) {}

  void draw() {
    *widget << CGAL::ORANGE;
    sdg.draw_skeleton(*widget);
  }

};

template< class T >
class Sites_layer : public CGAL::Qt_widget_layer {
private:
  T& sdg;

public:
  Sites_layer(T& sdg) : sdg(sdg) {}

  void draw() {
    *widget << CGAL::RED;
#if !defined (__POWERPC__)
    *widget << CGAL::PointSize(6);
    *widget << CGAL::LineWidth(3);
#endif
    {
      typename T::Finite_vertices_iterator vit;
      for (vit = sdg.finite_vertices_begin();
	   vit != sdg.finite_vertices_end(); ++vit) {
	typename T::Site_2 s = vit->site();
	*widget << CGAL::RED;
	if ( s.is_segment() ) {
	  *widget << s.segment();
	}
      }
    }
    {
      typename T::Finite_vertices_iterator vit;
      for (vit = sdg.finite_vertices_begin();
	   vit != sdg.finite_vertices_end(); ++vit) {
	typename T::Site_2 s = vit->site();
	if ( s.is_input() ) {
	  *widget << CGAL::RED;
	} else {
	  *widget << CGAL::YELLOW;
	}
	if ( s.is_point() ) {
	  *widget << s.point();
	}
      }
    }
  }
};


#endif // QT_LAYERS_H
