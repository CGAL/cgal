#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Delaunay_triangulation_2.h>


#include <qapplication.h>
#include <qmainwindow.h>
#include <qtoolbar.h>

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_widget_get_point.h>
#include <CGAL/IO/pixmaps/point.xpm>

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
    
    QToolBar  *tools_toolbar;
    tools_toolbar = new QToolBar("Tools", this, 
				QMainWindow::Top, TRUE, "Tools");
    addToolBar(tools_toolbar, Top, FALSE);
    get_point_but =  new QToolButton(QPixmap( (const char**)point_xpm ),
				  "Point Tool", 
				  0, 
				  this, 
				  SLOT(pointtool()), 
				  tools_toolbar, 
				  "Point Tool");
    get_point_but->setToggleButton(TRUE);
    widget->attach(&v);

    connect(widget, SIGNAL(new_cgal_object(CGAL::Object)), 
							this, SLOT(get_object(CGAL::Object)));
  }
  ~My_Window(){delete widget;}
private slots:
  //this function is called every time the toolbar button is pressed
  void pointtool(){
    if (get_point_but->isOn())
      widget->attach(&get_point);
    else
      widget->detach_current_tool();
  }

  //this function is called every time a tool creates a Cgal object
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
  CGAL::Qt_widget *widget;	//the instance of Qt_widget
  My_Layer v;		//an instance of a layer
  CGAL::Qt_widget_standard_toolbar *stoolbar; //the standard toolbar
  CGAL::Qt_widget_get_point<Rep> get_point;   
						//the generic tool that creates Cgal points
  QToolButton *get_point_but;	//the toolbar button
};

#include "tutorial8.moc"

int main( int argc, char **argv )
{
    QApplication app( argc, argv );
    My_Window W(600,600);
    app.setMainWidget( &W );
    W.show();
    W.setCaption("Using the Standard Toolbar");
    return app.exec();
}
