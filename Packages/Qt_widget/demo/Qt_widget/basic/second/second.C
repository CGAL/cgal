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
#include <CGAL/Delaunay_triangulation_2.h>


#include <qapplication.h>
#include <CGAL/IO/Qt_Widget.h>

typedef CGAL::Cartesian<double>		    Rep;
typedef CGAL::Point_2<Rep>		    Point;
typedef CGAL::Delaunay_triangulation_2<Rep> Delaunay;

Delaunay dt;

class My_Window : public CGAL::Qt_widget {
public:
  My_Window(int x, int y){ resize(x,y); };
private:
  
  //this event is called only when the user press mouse
  void mousePressEvent(QMouseEvent *e)
  {
    clear();
    dt.insert(Point(x_real(e->x()), y_real(e->y())));
    *this << dt;
  }
  
};

int main( int argc, char **argv )
{
    QApplication app( argc, argv );
    My_Window W(600,600);
    app.setMainWidget( &W );
    W.set_window(0, 600, 0, 600);
    W.show();
    return app.exec();
}

#endif // CGAL_USE_QT
