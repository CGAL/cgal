#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_Polygon_2.h>

//CGAL
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/min_quadrilateral_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Polygon_2_algorithms.h>

template <class R>
class Qt_layer_show_parallelogram : public CGAL::Qt_widget_layer
{
public:
  typedef typename R::Point_2	Point;
  typedef typename R::Segment_2	Segment;
  typedef typename R::Line_2	Line;


  typedef std::vector< Point >				    Container;
  typedef CGAL::Polygon_2< R, Container >		    Polygon;
  typedef CGAL::Creator_uniform_2< double, Point >	    Creator;

  Qt_layer_show_parallelogram(std::list<Point> *pl){list_of_points = pl;};

  void draw()
  {

    //Draw the MINIMUM PARALLELOGRAM
    widget->lock();
      Polygon pts;
      typedef typename std::list<Point>::const_iterator Listconstiter;
      for (Listconstiter i = (*list_of_points).begin(); i != (*list_of_points).end(); ++i)
        pts.push_back(Point(i->x(), i->y()));
      Polygon p;
      CGAL::convex_hull_points_2(
	pts.vertices_begin(), pts.vertices_end(), std::back_inserter(p));

      Polygon	kg;
      if (p.size() >= 3) {
	CGAL::min_parallelogram_2(
	  p.vertices_begin(), p.vertices_end(), std::back_inserter(kg));
      }
      RasterOp old = widget->rasterOp();	//save the initial raster mode
      widget->setRasterOp(XorROP);
      *widget << CGAL::FillColor(CGAL::GRAY);
      *widget << kg;
      widget->setRasterOp(old);
    widget->unlock();
  }
private:
  std::list<Point>  *list_of_points;  
};//end class 


template <class R>
class Qt_layer_show_strip : public CGAL::Qt_widget_layer
{
public:
  typedef typename R::Point_2	Point;
  typedef typename R::Segment_2	Segment;
  typedef typename R::Line_2	Line;


  typedef std::vector< Point >				    Container;
  typedef CGAL::Polygon_2< R, Container >		    Polygon;
  typedef CGAL::Creator_uniform_2< double, Point >	    Creator;

  Qt_layer_show_strip(std::list<Point> *pl){list_of_points = pl;};

  void draw()
  {
    //DRAW the Line Strips
    widget->lock();
      Polygon pts;
      typedef typename std::list<Point>::const_iterator Listconstiter;
      for (Listconstiter i = (*list_of_points).begin(); i != (*list_of_points).end(); ++i)
        pts.push_back(Point(i->x(), i->y()));
      Polygon p;
      CGAL::convex_hull_points_2(
	pts.vertices_begin(), pts.vertices_end(), std::back_inserter(p));

      std::vector< Line > ll;
      if (p.size() >= 2) 
	CGAL::min_strip_2(
	  p.vertices_begin(), p.vertices_end(), std::back_inserter(ll));
      *widget << CGAL::GREEN;
      for (typename std::vector<Line>::iterator it = ll.begin();
	   it != ll.end(); ++it)
	*widget << (*it);	    
    widget->unlock();

  }
private:
  std::list<Point>  *list_of_points;  
};//end class 


template <class R>
class Qt_layer_show_rectangle : public CGAL::Qt_widget_layer
{
public:
  typedef typename R::Point_2	Point;
  typedef typename R::Segment_2	Segment;
  typedef typename R::Line_2	Line;


  typedef std::vector< Point >				    Container;
  typedef CGAL::Polygon_2< R, Container >		    Polygon;
  typedef CGAL::Creator_uniform_2< double, Point >	    Creator;

  Qt_layer_show_rectangle(std::list<Point> *pl){list_of_points = pl;};

  void draw()
  {
    //Draw the MINIMUM RECTANGLE
    widget->lock();
      Polygon pts;
      typedef typename std::list<Point>::const_iterator Listconstiter;
      for (Listconstiter i = (*list_of_points).begin(); i != (*list_of_points).end(); ++i)
        pts.push_back(Point(i->x(), i->y()));
      Polygon p;
      CGAL::convex_hull_points_2(
	pts.vertices_begin(), pts.vertices_end(), std::back_inserter(p));
      
      Polygon kg;
      if (p.size() >= 3)
	CGAL::min_rectangle_2(
	  p.vertices_begin(), p.vertices_end(), std::back_inserter(kg));
      
      RasterOp old = widget->rasterOp();	//save the initial raster mode
      widget->setRasterOp(XorROP);
      *widget << CGAL::FillColor(CGAL::BLUE);      
      *widget << kg;
      widget->setRasterOp(old);  
    widget->unlock();

  }
private:
  std::list<Point>  *list_of_points;  
};//end class 



template <class R>
class Qt_layer_show_points : public CGAL::Qt_widget_layer
{
public:
  typedef typename R::Point_2	Point;
  
  Qt_layer_show_points(std::list<Point> *pl){list_of_points = pl;};

  void draw()
  {
    //Draw the points as CROSS
    widget->lock();
      *widget << CGAL::PointSize(3) << CGAL::PointStyle(CGAL::DISC);
      *widget << CGAL::GREEN;
      typename std::list<Point>::iterator itp = (*list_of_points).begin();
      while(itp!=(*list_of_points).end())
      {
	*widget << (*itp++);
      }
    widget->unlock();
  };//end draw();	
private:
  std::list<Point>  *list_of_points;  
};//end class 
