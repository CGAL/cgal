#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Delaunay_triangulation_2.h>


#include <qapplication.h>
#include <qmainwindow.h>

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>

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

class My_Widget : public CGAL::Qt_widget {
public:
  My_Widget(QMainWindow* c) : CGAL::Qt_widget(c) {};
private:
  //this event is called only when the user presses the mouse
  void mousePressEvent(QMouseEvent *e)
  {
    Qt_widget::mousePressEvent(e);
    dt.insert(Point(x_real(e->x()), y_real(e->y())));
    redraw();
  }
};

class My_Window : public QMainWindow{
public:
  My_Window(int x, int y)
  {
    widget = new My_Widget(this);
    setCentralWidget(widget);
    resize(x,y);
    widget->set_window(0, x, 0, y);
    
    //How to attach the standard toolbar
    stoolbar = new CGAL::Qt_widget_standard_toolbar(widget, this);
    this->addToolBar(stoolbar->toolbar(), Top, FALSE);
    setUsesBigPixmaps(true);

    widget->attach(&v);
  }
private:
  My_Widget *widget;
  My_Layer  v;
  CGAL::Qt_widget_standard_toolbar *stoolbar;
};

int main( int argc, char **argv )
{
    QApplication app( argc, argv );
    My_Window W(600,600);
    app.setMainWidget( &W );
    W.show();
    W.setCaption("Using the Standard Toolbar");
    return app.exec();
}

