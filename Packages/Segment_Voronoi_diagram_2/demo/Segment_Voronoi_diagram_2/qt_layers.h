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
    *widget << CGAL::BLUE;
    *widget << CGAL::PointSize(3);
    *widget << CGAL::LineWidth(3);
    svd.draw_primal(*widget);
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
    // svd.draw_skeleton(*widget);
    svd.draw_Voronoi_circles(*widget);
  }

};

template< class T >
class Sites_layer : public CGAL::Qt_widget_layer {
private:
  T& svd;

public:
  Sites_layer(T& svd) : svd(svd) {}

  void draw(){
    *widget << CGAL::RED;
    *widget << CGAL::PointSize(6);
    *widget << CGAL::LineWidth(3);
    svd.draw_sites(*widget);
  }
};


#endif // QT_LAYERS_H
