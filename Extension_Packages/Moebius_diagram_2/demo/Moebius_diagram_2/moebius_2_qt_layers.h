#ifndef CHR_MOEBIUS_2_QT_LAYERS_H
#define CHR_MOEBIUS_2_QT_LAYERS_H

#include <CGAL/Cartesian.h>
#include <CGAL/Conic_2.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_Conic_2.h>
#include <CGAL/IO/Qt_widget_Regular_triangulation_2.h>
#include <qobject.h>

template <class MD>
class Qt_layer_show_moebius : public CGAL::Qt_widget_layer
{
 public:
  Qt_layer_show_moebius (MD *md): _md (md) {}

  void draw () {
    if (! _md->constructed ()) return;
    *widget << CGAL::BLUE << *_md;
  }
 private:
  MD *_md;
};

template <class MD>
class Qt_layer_show_points : public CGAL::Qt_widget_layer
{
 public:
  typedef typename MD::Point Point;
  typedef typename MD::RT_3::Vertex_iterator Vertex_iterator;

  Qt_layer_show_points (MD *md): _md (md) {}
  void draw () {
    if (! _md->initialized ()) return;
    *widget << CGAL::GREEN << CGAL::PointSize (3) 
		<< CGAL::PointStyle (CGAL::DISC);    
    Vertex_iterator i = _md->rt().vertices_begin();
    Vertex_iterator stop = _md->rt().vertices_end();

    while (i != stop) {
      if (_md->is_finite (i))
	*widget << i->point ();
      ++i;
    }
  }
 private:
  MD *_md;
};

template <class MD>
class Qt_layer_show_lambdas : public CGAL::Qt_widget_layer
{
 public:
  typedef typename MD::Point Point;
  typedef typename MD::RT_3::Vertex_iterator Vertex_iterator;
  typedef typename MD::Circle Circle;

  Qt_layer_show_lambdas (MD *md): _md (md) {}
  void draw () {
    if (! _md->initialized ()) return;

    *widget << CGAL::RED;
    Vertex_iterator i = _md->rt().vertices_begin();
    Vertex_iterator stop = _md->rt().vertices_end();

    while (i != stop) {
      if (_md->is_finite (i))
	*widget << Circle (i->point (), 1 / (i->point().lambda()));
      ++i;
    }
  }
 private:
  MD *_md;
};


template <class MD>
class Qt_layer_show_mus : public CGAL::Qt_widget_layer
{
 public:
  typedef typename MD::Point Point;
  typedef typename MD::RT_3::Vertex_iterator Vertex_iterator;
  typedef typename MD::Circle Circle;

  Qt_layer_show_mus (MD *md): _md (md) {}
  void draw () {
    if (! _md->initialized ()) return;

    *widget << CGAL::GREEN;
    Vertex_iterator i = _md->rt().vertices_begin();
    Vertex_iterator stop = _md->rt().vertices_end();

    while (i != stop) {
      if (_md->is_finite (i))
	*widget << Circle (i->point (), i->point().mu());
      ++i;
    }
  }
 private:
  MD *_md;
};


#endif//CHR_MOEBIUS_2_QT_LAYERS_H
