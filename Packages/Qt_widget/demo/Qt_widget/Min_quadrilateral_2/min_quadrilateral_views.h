#include <CGAL/IO/Qt_widget_view.h>

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
class Qt_view_show_parallelogram : public CGAL::Qt_widget_view
{
public:
  typedef typename R::Point_2	Point;
  typedef typename R::Segment_2	Segment;
  typedef typename R::Line_2	Line;


  typedef CGAL::Polygon_traits_2< R>			    Traits;
  typedef std::vector< Point >				    Container;
  typedef CGAL::Polygon_2< Traits, Container >		    Polygon;
  typedef CGAL::Creator_uniform_2< double, Point >	    Creator;

  Qt_view_show_parallelogram(std::list<Point> *pl){list_of_points = pl;};

  void draw(CGAL::Qt_widget &win)
  {

    //Draw the MINIMUM PARALLELOGRAM
    win.lock();
      Polygon pts;
      typedef std::list<Point>::const_iterator Listconstiter;
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
      RasterOp old = win.rasterOp();	//save the initial raster mode
      win.setRasterOp(XorROP);
      win << CGAL::FillColor(CGAL::GRAY);
      win << kg;
      win.setRasterOp(old);
    win.unlock();
  }
private:
  std::list<Point>  *list_of_points;  
};//end class 


template <class R>
class Qt_view_show_strip : public CGAL::Qt_widget_view
{
public:
  typedef typename R::Point_2	Point;
  typedef typename R::Segment_2	Segment;
  typedef typename R::Line_2	Line;


  typedef CGAL::Polygon_traits_2< R>			    Traits;
  typedef std::vector< Point >				    Container;
  typedef CGAL::Polygon_2< Traits, Container >		    Polygon;
  typedef CGAL::Creator_uniform_2< double, Point >	    Creator;

  Qt_view_show_strip(std::list<Point> *pl){list_of_points = pl;};

  void draw(CGAL::Qt_widget &win)
  {
    //DRAW the Line Strips
    win.lock();
      Polygon pts;
      typedef std::list<Point>::const_iterator Listconstiter;
      for (Listconstiter i = (*list_of_points).begin(); i != (*list_of_points).end(); ++i)
        pts.push_back(Point(i->x(), i->y()));
      Polygon p;
      CGAL::convex_hull_points_2(
	pts.vertices_begin(), pts.vertices_end(), std::back_inserter(p));

      std::vector< Line > ll;
      if (p.size() >= 2) 
	CGAL::min_strip_2(
	  p.vertices_begin(), p.vertices_end(), std::back_inserter(ll));
      win << CGAL::GREEN;
      for (std::vector<Line>::iterator it = ll.begin(); it != ll.end(); ++it)
	win << (*it);	    
    win.unlock();

  }
private:
  std::list<Point>  *list_of_points;  
};//end class 


template <class R>
class Qt_view_show_rectangle : public CGAL::Qt_widget_view
{
public:
  typedef typename R::Point_2	Point;
  typedef typename R::Segment_2	Segment;
  typedef typename R::Line_2	Line;


  typedef CGAL::Polygon_traits_2< R>			    Traits;
  typedef std::vector< Point >				    Container;
  typedef CGAL::Polygon_2< Traits, Container >		    Polygon;
  typedef CGAL::Creator_uniform_2< double, Point >	    Creator;

  Qt_view_show_rectangle(std::list<Point> *pl){list_of_points = pl;};

  void draw(CGAL::Qt_widget &win)
  {
    //Draw the MINIMUM RECTANGLE
    win.lock();
      Polygon pts;
      typedef std::list<Point>::const_iterator Listconstiter;
      for (Listconstiter i = (*list_of_points).begin(); i != (*list_of_points).end(); ++i)
        pts.push_back(Point(i->x(), i->y()));
      Polygon p;
      CGAL::convex_hull_points_2(
	pts.vertices_begin(), pts.vertices_end(), std::back_inserter(p));
      
      Polygon kg;
      if (p.size() >= 3)
	CGAL::min_rectangle_2(
	  p.vertices_begin(), p.vertices_end(), std::back_inserter(kg));
      
      RasterOp old = win.rasterOp();	//save the initial raster mode
      win.setRasterOp(XorROP);
      win << CGAL::FillColor(CGAL::BLUE);      
      win << kg;
      win.setRasterOp(old);  
    win.unlock();

  }
private:
  std::list<Point>  *list_of_points;  
};//end class 



template <class R>
class Qt_view_show_points : public CGAL::Qt_widget_view
{
public:
  typedef typename R::Point_2	Point;
  
  Qt_view_show_points(std::list<Point> *pl){list_of_points = pl;};

  void draw(CGAL::Qt_widget &win)
  {
    
    
    //Draw the points as CROSS
    win.lock();
      win << CGAL::PointSize(7) << CGAL::PointStyle(CGAL::CROSS);
      win << CGAL::GREEN;
      std::list<Point>::iterator itp = (*list_of_points).begin();
      while(itp!=(*list_of_points).end())
      {
	win << (*itp++);
      }
    win.unlock();
  };//end draw();	
private:
  std::list<Point>  *list_of_points;  
};//end class 
