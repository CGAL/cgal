#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Delaunay_triangulation_2.h>


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

class My_Layer : public CGAL::Qt_widget_layer{
  void draw(CGAL::Qt_widget& widget){
    widget << CGAL::BLACK;
    widget << dt;
  }
};

class My_Window : public QMainWindow{
  Q_OBJECT
public:
  My_Window(int x, int y)
  {
    widget = new CGAL::Qt_widget(this);
    setCentralWidget(widget);
    resize(x,y);
    widget->show();
    widget->set_window(0, x, 0, y);
    
    //How to attach the standard toolbar
    stoolbar = new CGAL::Qt_widget_standard_toolbar(widget, this);
    this->addToolBar(stoolbar->toolbar(), Top, FALSE);
    
    widget->attach(&v);

    connect(widget, SIGNAL(new_cgal_object(CGAL::Object)), 
	  this, SLOT(get_object(CGAL::Object)));
    widget->attach(&get_point);
  }
  ~My_Window(){}
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
  My_Layer v;
  CGAL::Qt_widget_standard_toolbar *stoolbar;
  CGAL::Qt_widget_get_point<Rep> get_point;
};

#include "tutorial7.moc"

int main( int argc, char **argv )
{
    QApplication app( argc, argv );
    My_Window W(600,600);
    app.setMainWidget( &W );
    W.show();
    W.setCaption("Using the Standard Toolbar");
    return app.exec();
}
