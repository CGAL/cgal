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
#include <CGAL/IO/Qt_widget_standard_toolbar.h>

typedef CGAL::Cartesian<double>		    Rep;
typedef CGAL::Point_2<Rep>		    Point;
typedef CGAL::Delaunay_triangulation_2<Rep> Delaunay;

Delaunay dt;

class My_view : public CGAL::Qt_widget_view{
  void draw(CGAL::Qt_widget& win){
    win << CGAL::BLACK;
    win << dt;
  }
};

class My_tool : public CGAL::Qt_widget_tool{
public:
  My_tool(){};
private:
  void mousePressEvent(QMouseEvent *e)
  {
    if(e->button() == Qt::LeftButton)
    {
      double
	x=static_cast<double>(widget->x_real(e->x())),
	y=static_cast<double>(widget->y_real(e->y()));
      widget->new_object(CGAL::make_object(Point(x, y)));
    }
  }
};

class My_window : public QMainWindow{
  Q_OBJECT
public:
  My_window(int x, int y) : win(this)
  {
    setCentralWidget(&win);
    resize(x,y);
    win.show();
    win.set_window(0, x, 0, y);
    
    //How to attach the standard toolbar
    stoolbar = new CGAL::Standard_toolbar(&win, this);
    this->addToolBar(stoolbar->toolbar(), Top, FALSE);
    
    win.attach(&v);

    connect(&win, SIGNAL(new_cgal_object(CGAL::Object)), this, SLOT(get_object(CGAL::Object)));
    win.attach(t);
  }
private slots:
  void get_object(CGAL::Object obj)
  {
    Point p;
    if(CGAL::assign(p, obj))
    {
      dt.insert(p);
      win.redraw();
    }
  }
private:
  CGAL::Qt_widget win;
  My_view v;
  My_tool t;
  CGAL::Standard_toolbar *stoolbar;
};

#include "sixth.moc"

int main( int argc, char **argv )
{
    QApplication app( argc, argv );
    My_window W(600,600);
    app.setMainWidget( &W );
    W.show();
    W.setCaption("Using the Standard Toolbar");
    return app.exec();
}

#endif // CGAL_USE_QT
