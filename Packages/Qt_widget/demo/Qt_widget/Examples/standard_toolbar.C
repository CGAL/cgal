//demo/Qt_widget/Examples/Standard_toolbar.C
#ifndef CGAL_USE_QT
#include <iostream>
int main(int, char*){
  std::cout << "Sorry, this demo needs QT..." << std::endl; return 0;}
#else
#include <CGAL/Cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/IO/Qt_widget_Delaunay_triangulation_2.h>

#include <qapplication.h>
#include <qmainwindow.h>

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/point_generators_2.h>

typedef CGAL::Cartesian<double>             Rep;
typedef CGAL::Point_2<Rep>                  Point;
typedef CGAL::Delaunay_triangulation_2<Rep> Delaunay;

Delaunay dt;

class My_window : public QMainWindow{
  Q_OBJECT
public:
  My_window(int x, int y)
  {
    widget = new CGAL::Qt_widget(this);
    setCentralWidget(widget);
    resize(x,y);
    widget->show();
    widget->set_window(0, x, 0, y);

    CGAL::Random_points_in_disc_2<Point> g(500);
    for(int count=0; count<100; count++) {
      dt.insert(*g++);
    }
    
    //How to attach the standard toolbar
    std_toolbar = new CGAL::Qt_widget_standard_toolbar(widget, this,
						       "Standard Toolbar");

    connect(widget, SIGNAL(redraw_on_back()),
	    this, SLOT(redraw_win()) );
  }

private slots:	//functions
  void redraw_win()
  {
    *widget << dt;
  }

private:	//members
  CGAL::Qt_widget *widget;
  CGAL::Qt_widget_standard_toolbar *std_toolbar;
};

// moc_source_file: standard_toolbar.C
#include "standard_toolbar.moc"

int main( int argc, char **argv )
{
    QApplication app( argc, argv );
    My_window W(600,600);
    app.setMainWidget( &W );
    W.show();
    W.setCaption("Using the Standard Toolbar");
    return app.exec();
}
#endif
