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
#include <qmainwindow.h>

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_view.h>

typedef CGAL::Cartesian<double>		    Rep;
typedef CGAL::Point_2<Rep>		    Point;
typedef CGAL::Delaunay_triangulation_2<Rep> Delaunay;

Delaunay dt;

class My_view : public CGAL::Qt_widget_view{
  void draw(CGAL::Qt_widget& win){
    win << dt;
  }
};

class My_widget : public CGAL::Qt_widget {
public:
  My_widget(QMainWindow* c) : CGAL::Qt_widget(c) {};
private:
  //this event is called only when the user press mouse
  void mousePressEvent(QMouseEvent *e)
  {
    Qt_widget::mousePressEvent(e);
    dt.insert(Point(x_real(e->x()), y_real(e->y())));
    redraw();
  }
};

class My_window : public QMainWindow{
public:
  My_window(int x, int y) : win(this)
  {
    setCentralWidget(&win);
    resize(x,y);
    win.show();
    win.set_window(0, x, 0, y);
    win.attach(&v);
  }
private:
  My_widget win;
  My_view v;
};

int main( int argc, char **argv )
{
    QApplication app( argc, argv );
    My_window W(600,600);
    app.setMainWidget( &W );
    W.show();
    W.setCaption("Using QMainWindow QT class");
    return app.exec();
}

#endif // CGAL_USE_QT
