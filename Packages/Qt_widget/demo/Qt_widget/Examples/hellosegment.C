//demo/Qt_widget/Examples/hellosegment.C

#include <CGAL/Cartesian.h>
#include <CGAL/IO/Qt_widget.h>

#include <qapplication.h>

typedef CGAL::Cartesian<int> Rep;
typedef CGAL::Point_2<Rep> Point;
typedef CGAL::Segment_2<Rep> Segment;

int main( int argc, char **argv )
{
    QApplication app( argc, argv );
    CGAL::Qt_widget *W = new CGAL::Qt_widget();
    app.setMainWidget( W );
    W->resize(600, 600);
    W->set_window(0, 600, 0, 600);
    W->show();
    W->lock();
    *W << CGAL::BackgroundColor(CGAL::ORANGE) << CGAL::RED;
    *W << Segment(Point(100,100), Point(400,400));
    W->unlock();
    return app.exec();
}
