#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Delaunay_triangulation_2.h>


#include <qapplication.h>
#include <qmainwindow.h>

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/point_generators_2.h>

Delaunay dt;

class My_window : public QMainWindow{
public:
  My_window(int x, int y) : win(this)
  {
    setCentralWidget(&win);
    resize(x,y);
    win.show();
    win.set_window(0, x, 0, y);

    CGAL::Random_points_in_disc_2<Point> g(500);
    for(int count=0; count<100; count++) {
      dt.insert(*g++);
    }
    
    //How to attach the standard toolbar
    stoolbar = new CGAL::Standard_toolbar(&win, this);
    this->addToolBar(stoolbar->toolbar(), Top, FALSE);

    win.redraw();
  }
private:	//functions
  void redraw()
  {
    Qt_widget::redraw();
    *this << dt;
  }

private:	//members
  CGAL::Qt_widget win;
  CGAL::Standard_toolbar *stoolbar;
};

int main( int argc, char **argv )
{
    QApplication app( argc, argv );
    My_window W(600,600);
    app.setMainWidget( &W );
    W.show();
    W.setCaption("Using the Standard Toolbar");
    return app.exec();
}