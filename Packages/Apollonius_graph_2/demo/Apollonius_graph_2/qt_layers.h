#ifndef QT_LAYERS_H
#define QT_LAYERS_H

#include <CGAL/IO/Qt_widget_layer.h>

template< class T >
class Voronoi_diagram_layer : public CGAL::Qt_widget_layer {
private:
  T& ag;

public:
  Voronoi_diagram_layer(T& ag) : ag(ag) {}

  void draw() {
    *widget << CGAL::BLUE;
    ag.draw_dual(*widget);
  }

};

template< class T >
class Delaunay_graph_layer : public CGAL::Qt_widget_layer {
private:
  T& ag;

public:
  Delaunay_graph_layer(T& ag) : ag(ag) {}

  void draw(){
    *widget << CGAL::GREEN;
    ag.draw_primal(*widget);
  }
};

template< class T >
class Non_hidden_weighted_points_layer : public CGAL::Qt_widget_layer {
private:
  T& ag;

public:
  Non_hidden_weighted_points_layer(T& ag) : ag(ag) {}

  void draw(){
    *widget << CGAL::RED;
    ag.draw_non_hidden_weighted_points(*widget);
  }
};

template< class T >
class Hidden_weighted_points_layer : public CGAL::Qt_widget_layer {
private:
  T& ag;

public:
  Hidden_weighted_points_layer(T& ag) : ag(ag) {}

  void draw(){
    *widget << CGAL::Color(64,64,64);
    ag.draw_hidden_weighted_points(*widget);
  }
};

template< class T >
class Non_hidden_centers_layer : public CGAL::Qt_widget_layer {
private:
  T& ag;

public:
  Non_hidden_centers_layer(T& ag) : ag(ag) {}

  void draw(){
    *widget << CGAL::RED;
    ag.draw_non_hidden_weighted_point_centers(*widget);
  }
};

template< class T >
class Hidden_centers_layer : public CGAL::Qt_widget_layer {
private:
  T& ag;

public:
  Hidden_centers_layer(T& ag) : ag(ag) {}

  void draw(){
    *widget << CGAL::Color(64,64,64);
    ag.draw_hidden_weighted_point_centers(*widget);
  }
};

#endif // QT_LAYERS_H
