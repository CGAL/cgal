//demo/Qt_widget/Examples/tutorial2.C
#ifndef CGAL_USE_QT
#include <iostream>
int main(int, char*){
  std::cout << "Sorry, this demo needs QT..." << std::endl; return 0;}
#else
#include <CGAL/Cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/IO/Qt_widget_Triangulation_2.h>
#include <CGAL/IO/Qt_widget.h>
#include <qapplication.h>
#include <qmainwindow.h>

typedef CGAL::Cartesian<double>             K;
typedef K::Point_2                          Point;
typedef CGAL::Delaunay_triangulation_2<K>   Delaunay;

Delaunay dt;

class My_window : public QMainWindow {
  Q_OBJECT
public:
  My_window(int x, int y)
  {
    widget = new CGAL::Qt_widget(this);
    widget->resize(x,y);
    widget->set_window(0, x, 0, y);

    connect(widget, SIGNAL(redraw_on_back()),
	   this, SLOT(redraw_win()));

    connect(widget, SIGNAL(s_mousePressEvent(QMouseEvent*)),
	    this, SLOT(mousePressEvent(QMouseEvent*)));

    setCentralWidget(widget);
  };
private slots:  
  void redraw_win()
  {
    *widget << dt;
  }

  void mousePressEvent(QMouseEvent *e)
  {
    dt.insert(Point(widget->x_real(e->x()), widget->y_real(e->y())));
    widget->redraw();
  }

private: // private data member
  CGAL::Qt_widget* widget;
};

//moc_source_file : tutorial2.C
#include "tutorial2.moc"

int main( int argc, char **argv )
{
    QApplication app( argc, argv );
    My_window *w = new My_window(400,400);
    app.setMainWidget( w);
    w->show();
    return app.exec();
}
#endif
