//demo/Qt_widget/basic/tutorial2.C

#include <CGAL/Cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/IO/Qt_widget_Triangulation_2.h>
#include <CGAL/IO/Qt_widget.h>
#include <qapplication.h>

typedef CGAL::Cartesian<double>		    Rep;
typedef CGAL::Point_2<Rep>		    Point;
typedef CGAL::Delaunay_triangulation_2<Rep> Delaunay;

Delaunay dt;
typedef CGAL::Qt_widget Qt_widget;
class My_window : public Qt_widget {
  Q_OBJECT
public:
  My_window(int x, int y){
    resize(x,y);
    connect(this, SIGNAL(custom_redraw()),
	   this, SLOT(redraw_win()));
  };
private slots:  
  void redraw_win()
  {
    *this << dt;
  }
private:
  //this event is called only when the user presses the mouse
  void mousePressEvent(QMouseEvent *e)
  {
    Qt_widget::mousePressEvent(e);
    dt.insert(Point(x_real(e->x()), y_real(e->y())));
    redraw();
  }
};

//moc_source_file : tutorial2.C
#include "tutorial2.moc"

int main( int argc, char **argv )
{
    QApplication app( argc, argv );
    My_window *W = new My_window(600,600);
    app.setMainWidget( W );
    W->show();
    W->set_window(0, 600, 0, 600);
    return app.exec();
}
