//demo/Qt_widget/Examples/layer.C

#include <CGAL/Cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/IO/Qt_widget_Delaunay_triangulation_2.h>
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_get_point.h>

#include <qapplication.h>

typedef CGAL::Cartesian<double>             Rep;
typedef CGAL::Point_2<Rep>                  Point;
typedef CGAL::Delaunay_triangulation_2<Rep> Delaunay;
typedef CGAL::Qt_widget						Qt_widget;

Delaunay dt;

class My_Layer : public CGAL::Qt_widget_layer{
  void draw(CGAL::Qt_widget& win){
    win << dt;
  }
};

class My_Window : public Qt_widget {
  Q_OBJECT
public:
  My_Window(int x, int y){
    resize(x,y);
    attach(&get_point);
    attach(&v);
    connect(this, SIGNAL(new_cgal_object(CGAL::Object)), 
            this, SLOT(get_new_object(CGAL::Object)));
  };
private:	//members
  CGAL::Qt_widget_get_point<Rep> get_point;
  My_Layer v;
private slots:
  void get_new_object(CGAL::Object obj)
  {
    Point p;
    if (CGAL::assign(p, obj)) { 
      dt.insert(p);
    }
    redraw();
  }

}; //endclass

//  moc_source_file : layer.C
#include "layer.moc"

int main( int argc, char **argv )
{
    QApplication app( argc, argv );
    My_Window *W = new My_Window(600,600);
    app.setMainWidget( W );
    W->set_window(0, 600, 0, 600);
    W->show();
    return app.exec();
}
