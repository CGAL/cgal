#ifndef QT_LAYERS_H
#define QT_LAYERS_H

#include <CGAL/IO/Qt_widget_layer.h>
#include "Virtual_Voronoi_diagram_2.h"

template< class T >
class Voronoi_diagram_layer : public CGAL::Qt_widget_layer
{
private:
  T* vd;

public:
  Voronoi_diagram_layer(T* vd) : vd(vd) {}

  void set(T* vvd) { vd = vvd; }

  void draw() {
    *widget << CGAL::BLUE;
#if !defined (__POWERPC__)
    *widget << CGAL::PointSize(3);
    *widget << CGAL::LineWidth(3);
#endif
    vd->draw_diagram(*widget);
  }
};

template< class T >
class Sites_layer : public CGAL::Qt_widget_layer
{
private:
  T* vd;

public:
  Sites_layer(T* vd) : vd(vd) {}

  void set(T* vvd) { vd = vvd; }

  void draw() {
    *widget << CGAL::RED;
#if !defined (__POWERPC__)
    *widget << CGAL::PointSize(6);
    *widget << CGAL::LineWidth(3);
#endif
    vd->draw_sites(*widget);
  }
};


#endif // QT_LAYERS_H
