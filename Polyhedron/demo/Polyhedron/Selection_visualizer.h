#ifndef SELECTION_VISUALIZER_H
#define SELECTION_VISUALIZER_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Viewer_interface.h>
#include <QPainter>

// Class for visualizing selection
// provides mouse selection functionality
class Q_DECL_EXPORT Selection_visualizer
{

 private:
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef K::Point_2 Point_2;
  typedef K::Point_3 Point_3;
  typedef CGAL::Polygon_2<K> Polygon_2;
  typedef std::vector<Point_2> Polyline_2;
  typedef std::vector<Polyline_2> Polylines;
  typedef CGAL::Three::Scene_item::Bbox Bbox;

  bool rectangle;
  std::vector<Point_2> contour_2d;
  Polylines* polyline;
  Bbox point_set_bbox;
  Polygon_2 domain_freeform;

public:
  CGAL::Bbox_2 domain_rectangle;

  Selection_visualizer(bool rectangle, const Bbox& point_set_bbox)
    : rectangle (rectangle), point_set_bbox (point_set_bbox)
  {
    polyline = new Polylines(0);
    polyline->push_back(Polyline_2());
  }
  ~Selection_visualizer() {
  }

  void render(QImage& image) const {

    CGAL::Three::Viewer_interface* viewer = static_cast<CGAL::Three::Viewer_interface*>(*CGAL::QGLViewer::QGLViewerPool().begin());

    QPen pen;
    pen.setColor(QColor(Qt::green));
    pen.setWidth(5);
    QImage temp(image);

    QPainter *painter = new QPainter(&temp);
    painter->setPen(pen);
    for(std::size_t i=0; i<polyline->size(); ++i)
    {
      Polyline_2 poly = (*polyline)[i];
      if(!poly.empty())
        for(std::size_t j=0; j<poly.size()-1; ++j)
        {
          painter->drawLine(poly[j].x(), poly[j].y(), poly[j+1].x(), poly[j+1].y());
        }
    }
    painter->end();
    delete painter;
    viewer->set2DSelectionMode(true);
    viewer->setStaticImage(temp);
    viewer->update();
  }

  Polyline_2& poly() const
  { return polyline->front(); }

  bool update_polyline () const
  {
    if (contour_2d.size() < 2 ||
        (!(poly().empty()) && contour_2d.back () == poly().back()))
      return false;

    if (rectangle)
      {
        poly().clear();

        poly().push_back ( Point_2 (domain_rectangle.xmin(),
                                domain_rectangle.ymin()));
        poly().push_back ( Point_2 (domain_rectangle.xmax(),
                                domain_rectangle.ymin()));
        poly().push_back ( Point_2 (domain_rectangle.xmax(),
                                domain_rectangle.ymax()));
        poly().push_back ( Point_2 (domain_rectangle.xmin(),
                                domain_rectangle.ymax()));
        poly().push_back ( Point_2 (domain_rectangle.xmin(),
                                                domain_rectangle.ymin()));

      }
    else
      {
        if (!(poly().empty()) && contour_2d.back () == poly().back())
          return false;

        poly().clear();

        for (unsigned int i = 0; i < contour_2d.size (); ++ i)
          poly().push_back (contour_2d[i]);
      }
    return true;
  }


  void sample_mouse_path(QImage& image)
  {
    CGAL::QGLViewer* viewer = *CGAL::QGLViewer::QGLViewerPool().begin();
    const QPoint& p = viewer->mapFromGlobal(QCursor::pos());

    if (rectangle && contour_2d.size () == 2)
      {
        contour_2d[1] = Point_2 (p.x (), p.y ());
        domain_rectangle = CGAL::bbox_2 (contour_2d.begin (), contour_2d.end ());
      }
    else
      contour_2d.push_back (Point_2 (p.x (), p.y ()));

    if (update_polyline ())
      {

        render(image);
      }
  }

  void apply_path()
  {
    update_polyline ();
    domain_rectangle = CGAL::bbox_2 (contour_2d.begin (), contour_2d.end ());
    if (!rectangle)
      domain_freeform = Polygon_2 (contour_2d.begin (), contour_2d.end ());
  }

  bool is_selected (CGAL::qglviewer::Vec& p)
  {
    if (domain_rectangle.xmin () < p.x &&
        p.x < domain_rectangle.xmax () &&
        domain_rectangle.ymin () < p.y &&
        p.y < domain_rectangle.ymax ())
      {
        if (rectangle)
          return true;
/*
 * domain_freeform.has_on_bounded_side() requires the polygon to be simple, which is never the case.
 * However, it works very well even if the polygon is not simple, so we use this instead to avoid
 * the cgal_assertion on is_simple().*/


        if (CGAL::bounded_side_2(domain_freeform.container().begin(),
                                 domain_freeform.container().end(),
                                 Point_2(p.x, p.y),
                                 domain_freeform.traits_member())  == CGAL::ON_BOUNDED_SIDE)
          return true;
      }
    return false;
  }


}; // end class Selection_visualizer
#endif // SELECTION_VISUALIZER_H
