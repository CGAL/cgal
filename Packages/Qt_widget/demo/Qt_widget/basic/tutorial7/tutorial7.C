//demo/Qt_widget/basic/tutorial7/tutorial7.C

#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/IO/Qt_widget_Delaunay_triangulation_2.h>

#include <qapplication.h>
#include <qmainwindow.h>

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_widget_get_point.h>

typedef CGAL::Cartesian<double>		    Rep;
typedef CGAL::Point_2<Rep>		    Point;
typedef CGAL::Delaunay_triangulation_2<Rep> Delaunay;

Delaunay dt;

class My_layer : public CGAL::Qt_widget_layer{
  void draw(){
    *widget << CGAL::BLACK;
    *widget << dt;
  }
};

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
    
    //How to attach the standard toolbar
    std_toolbar = new CGAL::Qt_widget_standard_toolbar(widget, this);
    this->addToolBar(std_toolbar->toolbar(), Top, FALSE);
    
    widget->attach(&v);

    connect(widget, SIGNAL(new_cgal_object(CGAL::Object)), 
	  this, SLOT(get_object(CGAL::Object)));
    widget->attach(&get_point);
  }
  ~My_window(){}
private slots:
  void get_object(CGAL::Object obj)
  {
    Point p;
    if(CGAL::assign(p, obj))
    {
      dt.insert(p);
      widget->redraw();
      *widget << CGAL::RED << p;
    }
  }
private:
  CGAL::Qt_widget *widget;
  My_layer v;
  CGAL::Qt_widget_standard_toolbar *std_toolbar;
  CGAL::Qt_widget_get_point<Rep> get_point;
};

#include "tutorial7.moc"

int main( int argc, char **argv )
{
    QApplication app( argc, argv );
    My_window W(600,600);
    app.setMainWidget( &W );
    W.show();
    W.setCaption("Using the Standard Toolbar");
    return app.exec();
}
