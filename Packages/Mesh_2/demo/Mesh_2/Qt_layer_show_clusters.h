// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// file          : 
// package       : 
// author(s)     : Laurent Rineau
// release       : 
// release_date  : 
//
//
//
// ============================================================================

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>

class Show_clusters_aux : public CGAL::Qt_widget_layer
{
  Q_OBJECT
private:
  virtual voir reinit_clusters() {}

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
  typedef typename Conform::Triangulation       DT;
  typedef typename Conform::Segment		Segment;
  typedef typename Conform::Face_handle		Face_handle;
  typedef typename Conform::Vertex_handle	Vertex_handle;
  typedef typename Conform::Geom_traits::FT	FT;
  typedef typename Conform::Cluster_vertices_iterator CVIt;
  typedef typename Conform::Vertices_in_cluster_iterator ViCIt
  typedef std::list<Point> List_of_points;
  typedef typename List_of_points::const_iterator Point_iterator;

  Show_clusters(Conform &conform,
			 Color c = CGAL::GREEN,
			 int pointsize = 3,
			 PointStyle pointstyle = CGAL::DISC,
			 Color lc = CGAL::RED,
			 int linewidth = 2)
    : c(conform), first_time(true),  dt(), _color(c),
      size(pointsize), style(pointstyle), _line_color(lc),
      width(linewidth) {}

  void reinit_clusters()
  {
    dt.clear();

    for(CVIt it = c.clusters_vertices_begin();
	it != c.clusters_vertices_end();
	++it)
      dt.push_back( (*it)->point() );
  }

  void draw(){first_time = true;}

  void mouseMoveEvent(QMouseEvent *e)
  {
    if (dt.dimension()<1) return;

    FT x, y;
    widget->x_real(e->x(), x);
    widget->y_real(e->y(), y);
    Point p(x, y);

    RasterOp old = widget->rasterOp();	//save the initial raster mode
    widget->setRasterOp(XorROP);
    widget->lock();

    Vertex_handle v = dt.nearest_vertex(p);
    *widget << _color << CGAL::PointSize (pointsize)
	    << CGAL::PointStyle (pointstyle);

    if(!first_time)
      *widget << oldPoint;
    *widget << v->point();

    *widget << _line_color << CGAL::LineWidth(width);
    if(!first_time)
      for(Point_iterator pIt = oldPoints.begin();
	  pIt != oldPoints.end();
	  ++pIt)
	*widget << Segment(oldPoint, *pIt);
    oldPoints.clear();

    typename Conform::Locate_type lt;
    int i;
    Face_handle fh = c.locate(v->point(), lt, t);
    CGAL_assertion( lt == this->VERTEX );

    Vertex_handle v2 = fh->vertex(i);

    int n = c.number_of_clusters_at_vertex(v2);

    for(int j = 0; j < n; ++j)
      {
	std::pair<ViCIt> seq = c.vertices_in_cluster_sequence(v2, j);
	for(ViCIt it = seq.first;
	    it != seq.second;
	    ++it)
	  {
	    oldPoints.push_back((*it)->point());
	    *widget << Segment(*v2, (*it)->point());
	  }
      }

    widget->unlock();
    widget->setRasterOp(old);
    oldPoint = v->point();
    first_time = false;
  }

  void leaveEvent(QEvent *)
  {
    widget->lock();
    RasterOp old = widget->rasterOp();	//save the initial raster mode
    widget->setRasterOp(XorROP);

    *widget << _color << CGAL::PointSize (pointsize) 
	    << CGAL::PointStyle (pointstyle);
    if(!first_time) *widget << oldPoint;

    widget->unlock();
    widget->setRasterOp(old);
    first_time = true;
  }

private:
  Conform& c;
  DT dt;
  Point oldPoint, newPoint;
  List_of_points oldPoints;
  bool  first_time;
  Color _color;
  int size;
  PointStyle style;
  Color _line_color;
  int width;
};

// moc_source_file: Show_clusters.h
#include "Show_clusters.moc"
