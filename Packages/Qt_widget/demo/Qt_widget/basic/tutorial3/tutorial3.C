#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <qapplication.h>

typedef CGAL::Cartesian<double>		    Rep;
typedef CGAL::Point_2<Rep>		    Point;
typedef CGAL::Delaunay_triangulation_2<Rep> Delaunay;

Delaunay dt;

class My_Layer : public CGAL::Qt_widget_layer{
  void draw(CGAL::Qt_widget& win){
    win << dt;
  }
};

class My_Window : public CGAL::Qt_widget {
public:
  My_Window(int x, int y)
  {
    resize(x,y);
    attach(&v);
  };
private:
  //this event is called only when the user presses the mouse
  void mousePressEvent(QMouseEvent *e)
  {
    Qt_widget::mousePressEvent(e);
    dt.insert(Point(x_real(e->x()), y_real(e->y())));
    redraw();
  }
  My_Layer v;
};

int main( int argc, char **argv )
{
    QApplication app( argc, argv );
    My_Window *W = new My_Window(600,600);
    app.setMainWidget(W);
    W->show();
    W->set_window(0, 600, 0, 600);
    return app.exec();
}
