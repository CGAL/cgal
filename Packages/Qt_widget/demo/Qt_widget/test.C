// if QT is not installed, a message will be issued in runtime.
#ifndef CGAL_USE_QT
#include <iostream>

int main(int, char*)
{

  std::cout << "Sorry, this demo needs QT...";
  std::cout << std::endl;

  return 0;
}

#else

#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Bbox_2.h>
#include <list>
#include <CGAL/Polygon_2.h>
#include <qapplication.h>
#include <CGAL/IO/Qt_Window.h>

typedef CGAL::Cartesian<int> Rep;
typedef CGAL::Point_2<Rep> Point;
typedef CGAL::Circle_2<Rep> Circle;
typedef CGAL::Segment_2<Rep> Segment;
typedef CGAL::Line_2<Rep> Line;
typedef CGAL::Ray_2<Rep> Ray;
typedef CGAL::Triangle_2<Rep> Triangle;
typedef CGAL::Iso_rectangle_2<Rep> Rectangle;
typedef CGAL::Bbox_2 BBox;
typedef std::list<Point> Container;
typedef CGAL::Polygon_2<Rep,Container> Polygon;

class My_Window : public CGAL::Qt_widget {
public:
  My_Window(int x, int y) { resize(x,y); };
private:
  void	mousePressEvent(QMouseEvent *e)
  {
    using namespace CGAL;
    if(e->button() == Qt::RightButton)
      qApp->quit();
    else
      {
    this->lock();
    *this << BackgroundColor(ORANGE) << RED <<
	  LineWidth(3) << PointSize(3) << PointStyle(DISC);
    *this << Segment(Point(10,20),Point(300,400));
    *this << LineWidth(5) << GREEN << FillColor(BLACK) <<
      Circle(Point(400,400),50*50);
    *this << LineWidth(1) << noFill << Circle(Point(300,300),300*300);
    *this << BLUE << LineWidth(2);
    *this << Segment(Point(200,200),Point(400,400));
    *this << Segment(Point(200,400),Point(400,200));
    this->setFilled(true);
    *this << RED << Triangle(Point(150,300),
				   Point(150,350),
				   Point(100,325));
    *this << FillColor(RED) << Rectangle(Point(320,220),
					       Point(350,240));
    *this << DEEPBLUE << BBox(100,80,260,140);
    Polygon p;
    p.push_back(Point(300,30));
    p.push_back(Point(400,30));
    p.push_back(Point(500,130));
    p.push_back(Point(400,180));
    p.push_back(Point(300,130));
    *this << p;
    *this << PointSize(2) << Point(10,10);
    *this << PointSize(10)
	  << PointStyle(PIXEL) << Point(120,100)
	  << PointStyle(CROSS) << Point(140,100)
	  << PointStyle(PLUS) << Point(160,100)
	  << PointStyle(CIRCLE) << Point(180,100)
	  << PointStyle(DISC) << Point(200,100)
	  << PointStyle(RECT) << Point(220,100)
	  << PointStyle(BOX) << Point(240,100);
    *this << PointSize(5)
	  << PointStyle(PIXEL) << Point(120,120)
	  << PointStyle(CROSS) << Point(140,120)
	  << PointStyle(PLUS) << Point(160,120)
	  << PointStyle(CIRCLE) << Point(180,120)
	  << PointStyle(DISC) << Point(200,120)
	  << PointStyle(RECT) << Point(220,120)
	  << PointStyle(BOX) << Point(240,120);
    *this << LineWidth(1) << RED <<
      Line(Point(310,320), Point(320,340)) <<
      Line(Point(300,300), Point(150,350));
    *this << Ray(Point(200,400), Point(180,430))
	  << Ray(Point(200,400), Point(180,370));
    this->unlock();
      }
  };
};

int main( int argc, char **argv )
{
    QApplication app( argc, argv );
    My_Window   W(600,600);
    app.setMainWidget( &W );
    W.show();
    W << Segment(Point(250,200),Point(350,400));
    W << Segment(Point(200,350),Point(400,250));
    return app.exec();
}

#endif // CGAL_USE_QT
