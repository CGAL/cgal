//demo/Qt_widget/Examples/hellosegment.C
#ifndef CGAL_USE_QT
#include <iostream>
int main(int, char*){
  std::cout << "Sorry, this demo needs QT..." << std::endl; return 0;}
#else
#include <CGAL/Cartesian.h>
#include <CGAL/IO/Qt_widget.h>

#include <qapplication.h>

typedef CGAL::Cartesian<int> Rep;
typedef CGAL::Point_2<Rep> Point_2;
typedef CGAL::Segment_2<Rep> Segment;

int main( int argc, char **argv )
{
  QApplication app( argc, argv );
  CGAL::Qt_widget *w = new CGAL::Qt_widget();
  app.setMainWidget( w );
  w->resize(600, 600);
  w->set_window(0, 600, 0, 600);
  w->show();
  w->lock();
  *w << CGAL::BackgroundColor(CGAL::ORANGE) << CGAL::RED;
  *w << Segment(Point_2(100,100), Point_2(400,400));
  w->unlock();
  return app.exec();
}
#endif
