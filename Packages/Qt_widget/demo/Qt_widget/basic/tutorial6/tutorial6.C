//demo/Qt_widget/basic/tutorial6/tutorial6.C

#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/IO/Qt_widget_Delaunay_triangulation_2.h>

#include <qapplication.h>
#include <qmainwindow.h>

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>

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

class My_input_layer : public CGAL::Qt_widget_layer{
public:
  My_input_layer(){};
private:
  void mousePressEvent(QMouseEvent *e)
  {
    if(e->button() == Qt::LeftButton)
    {
      double x=static_cast<double>(widget->x_real(e->x()));
      double y=static_cast<double>(widget->y_real(e->y()));
      widget->new_object(CGAL::make_object(Point(x, y)));
    }
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
    widget->set_window(0, x, 0, y);
    
    //How to attach the standard toolbar
    std_toolbar = new CGAL::Qt_widget_standard_toolbar(widget, this);
    this->addToolBar(std_toolbar->toolbar(), Top, FALSE);
    
    widget->attach(&v);

    connect(widget, SIGNAL(new_cgal_object(CGAL::Object)), 
	    this, SLOT(get_object(CGAL::Object)));
    widget->attach(&t);
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
    }
  }
private:
  CGAL::Qt_widget *widget;
  My_layer v;
  My_input_layer t;
  CGAL::Qt_widget_standard_toolbar *std_toolbar;
};

#include "tutorial6.moc"

int main( int argc, char **argv )
{
    QApplication app( argc, argv );
    My_window W(600,600);
    app.setMainWidget( &W );
    W.show();
    W.setCaption("Using the Standard Toolbar");
    return app.exec();
}
