//demo/Qt_widget/basic/tutorial1/tutorial1.C

#include <CGAL/Cartesian.h>
#include <CGAL/Bbox_2.h>
#include <list>
#include <CGAL/Polygon_2.h>
#include <CGAL/IO/Qt_widget_Polygon_2.h>
#include <qapplication.h>
#include <CGAL/IO/Qt_widget.h>

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

int main( int argc, char **argv )
{
    QApplication app( argc, argv );
    using namespace CGAL;
    CGAL::Qt_widget * W = new CGAL::Qt_widget();
    app.setMainWidget( W );
    W->resize(600, 600);
    W->set_window(0, 600, 0, 600);
    W->show();
    //painting something on the screen
    W->lock();
    
    *W << BackgroundColor(ORANGE) << RED <<
	  LineWidth(3) << PointSize(3) << PointStyle(DISC);
    *W << Segment(Point(10,20),Point(300,400));
    *W << LineWidth(5) << GREEN << FillColor(BLACK) <<
      Circle(Point(400,400),50*50);
    *W << LineWidth(1) << noFill << Circle(Point(300,300),300*300);
    *W << BLUE << LineWidth(2);
    *W << Segment(Point(200,200),Point(400,400));
    *W << Segment(Point(200,400),Point(400,200));
    W->setFilled(TRUE);
    *W << RED << Triangle(Point(150,300),
				   Point(150,350),
				   Point(100,325));
    *W << FillColor(RED) << Rectangle(Point(320,220),
					       Point(350,240));
    *W << DEEPBLUE << BBox(100,80,260,140);
    Polygon p;
    p.push_back(Point(300,30));
    p.push_back(Point(400,30));
    p.push_back(Point(500,130));
    p.push_back(Point(400,180));
    p.push_back(Point(300,130));
    *W << p;
    *W << Ray(Point(200,400), Point(180,430))
      << Ray(Point(200,400), Point(180,370));
    W->unlock();

    return app.exec();
}
