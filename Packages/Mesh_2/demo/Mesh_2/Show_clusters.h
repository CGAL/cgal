#ifndef MESH_2_SHOW_CLUSTERS_H
#define MESH_2_SHOW_CLUSTERS_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>

#include <CGAL/Delaunay_triangulation_2.h>

#include <qpixmap.h>

class Show_clusters_aux : public CGAL::Qt_widget_layer
{
  Q_OBJECT
private:
  virtual void reinit_clusters() {}

public slots:
  void reinitClusters()
  {
    reinit_clusters();
  }
};

template <class Conform>
class Show_clusters : public Show_clusters_aux
{
public:
  typedef typename Conform::Point		Point;
  typedef typename Conform::Geom_traits         Geom_traits;
  typedef CGAL::Delaunay_triangulation_2<Geom_traits> DT;
  typedef typename DT::Finite_vertices_iterator Vertices_iterator;
  typedef typename DT::Vertex_handle DT_vertex_handle;
  typedef typename Conform::Segment		Segment;
  typedef typename Conform::Face_handle		Face_handle;
  typedef typename Conform::Vertex_handle	Vertex_handle;
  typedef typename Conform::Geom_traits::FT	FT;
  typedef typename Conform::Cluster_vertices_iterator CVIt;
  typedef typename Conform::Vertices_in_cluster_iterator ViCIt;
  typedef std::list<Point> List_of_points;
  typedef typename List_of_points::const_iterator Point_iterator;

  Show_clusters(Conform &conform,
		CGAL::Color color = CGAL::GREEN,
		int pointsize = 3,
		CGAL::PointStyle pointstyle = CGAL::DISC,
		CGAL::Color lc = CGAL::RED,
		int linewidth = 2)
    : c(conform), dt(), _color(color),
      size(pointsize), style(pointstyle), _line_color(lc),
      width(linewidth) 
    {
      for(CVIt it = conform.clusters_vertices_begin();
	  it != conform.clusters_vertices_end();
	  ++it)
	dt.push_back( (*it)->point() );
    }

  void activating()
  {
    reinit_clusters();
  }

  void reinit_clusters()
  {
    if(!is_active()) return;

    dt.clear();

    for(CVIt it = c.clusters_vertices_begin();
	it != c.clusters_vertices_end();
	++it)
      {
	dt.push_back( (*it)->point() );
      }
  }

  void draw()
  {
    widget->lock();

    QColor oldColor = widget->color();
    int oldPointSize = widget->pointSize();
    CGAL::PointStyle oldStyle = widget->pointStyle();

    *widget << _color << CGAL::PointStyle(style)
	    << CGAL::PointSize(size);

    for(Vertices_iterator it = dt.finite_vertices_begin();
	it != dt.finite_vertices_end();
	++it)
      *widget << it->point();

    widget->setPointStyle(oldStyle);
    widget->setPointSize(oldPointSize);
    widget->setColor(oldColor);

    widget->unlock();
    oldPixmap = widget->get_pixmap();
  }

  void mouseMoveEvent(QMouseEvent *e)
  {
    FT x, y;
    widget->x_real(e->x(), x);
    widget->y_real(e->y(), y);
    Point p(x, y);

    DT_vertex_handle v = dt.nearest_vertex(p);

    if(v == NULL) return;
    if(v == oldVertex) return;

    oldVertex = v;

    QColor oldColor = widget->color();
    int oldWidth = widget->lineWidth();
    
    widget->lock();

    widget->get_painter().drawPixmap(0, 0, oldPixmap);

    *widget << _line_color << CGAL::LineWidth(width);

    typename Conform::Locate_type lt;
    int i;
    Face_handle fh = c.locate(v->point(), lt, i);
    CGAL_assertion( lt == Conform::VERTEX );

    Vertex_handle v2 = fh->vertex(i);

    int n = c.number_of_clusters_at_vertex(v2);

    for(int j = 0; j < n; ++j)
      {
	std::pair<ViCIt,ViCIt> seq = c.vertices_in_cluster_sequence(v2, j);
	for(ViCIt it = seq.first;
	    it != seq.second;
	    ++it)
	  *widget << Segment(v2->point(), (*it)->point());
      }

    widget->setLineWidth(oldWidth);
    widget->setColor(oldColor);

    widget->unlock();
  }

  void leaveEvent(QEvent *)
  {
    widget->get_painter().drawPixmap(0, 0, oldPixmap);
    widget->update();
  }

private:
  Conform& c;
  DT dt;
  DT_vertex_handle oldVertex;
  QPixmap oldPixmap;
  bool  should_restore_pixmap;
  CGAL::Color _color;
  int size;
  CGAL::PointStyle style;
  CGAL::Color _line_color;
  int width;
};

#endif
