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
class Visible_sites_layer : public CGAL::Qt_widget_layer {
private:
  T& ag;

public:
  Visible_sites_layer(T& ag) : ag(ag) {}

  void draw(){
    *widget << CGAL::RED;

    typename T::Visible_sites_iterator it;
    for (it = ag.visible_sites_begin();
	 it != ag.visible_sites_end(); it++) {
      //      *widget << to_circle(*it);
      //      *widget << it->point();
      *widget << *it;
    }
  }
};

template< class T >
class Hidden_sites_layer : public CGAL::Qt_widget_layer {
private:
  T& ag;

public:
  Hidden_sites_layer(T& ag) : ag(ag) {}

  void draw(){
    *widget << CGAL::Color(64,64,64);

    typename T::Hidden_sites_iterator it;
    for (it = ag.hidden_sites_begin();
	 it != ag.hidden_sites_end(); it++) {
      //      *widget << to_circle(*it);
      //      *widget << it->point();
      *widget << *it;
    }
  }
};


#endif // QT_LAYERS_H
