#ifndef QT_LAYERS_H
#define QT_LAYERS_H

#include <CGAL/IO/Qt_widget_layer.h>

template< class T >
class Voronoi_diagram_layer : public CGAL::Qt_widget_layer {
private:
  T& svd;

public:
  Voronoi_diagram_layer(T& svd) : svd(svd) {}

  void draw() {
#if 1
    draw_diagram(*widget, svd);
#else 

    *widget << CGAL::BLUE;
#if !defined (__POWERPC__)
    *widget << CGAL::PointSize(3);
    *widget << CGAL::LineWidth(3);
#endif
    svd.draw_dual(*widget);

#endif
  }

};

template< class T >
class Skeleton_layer : public CGAL::Qt_widget_layer {
private:
  T& svd;

public:
  Skeleton_layer(T& svd) : svd(svd) {}

  void draw() {
    *widget << CGAL::ORANGE;
    svd.draw_skeleton(*widget);
  }

};

template< class T >
class Sites_layer : public CGAL::Qt_widget_layer {
private:
  T& svd;

public:
  Sites_layer(T& svd) : svd(svd) {}

  void draw() {
    *widget << CGAL::RED;
#if !defined (__POWERPC__)
    *widget << CGAL::PointSize(6);
    *widget << CGAL::LineWidth(3);
#endif
    {
      typename T::Finite_vertices_iterator vit;
      for (vit = svd.finite_vertices_begin();
	   vit != svd.finite_vertices_end(); ++vit) {
	typename T::Site_2 s = vit->site();
	*widget << CGAL::RED;
	if ( s.is_segment() ) {
	  *widget << s.segment();
	}
      }
    }
    {
      typename T::Finite_vertices_iterator vit;
      for (vit = svd.finite_vertices_begin();
	   vit != svd.finite_vertices_end(); ++vit) {
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
