//demo/Qt_widget/basic/tutorial8/tutorial8.C

#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/IO/Qt_widget_Delaunay_triangulation_2.h>

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
    
    QToolBar  *layers_toolbar;
    layers_toolbar = new QToolBar("Tools", this, 
				QMainWindow::Top, TRUE, "Tools");
    addToolBar(layers_toolbar, Top, FALSE);
    get_point_button = new QToolButton(layers_toolbar, "Get Point");
    get_point_button->setPixmap(QPixmap( (const char**)point_xpm ));
    get_point_button->setToggleButton(TRUE);
    widget->attach(&v);
    widget->attach(&get_point);

    connect(get_point_button, SIGNAL(stateChanged(int)),
	    &get_point, SLOT(stateChanged(int)));

    connect(widget, SIGNAL(new_cgal_object(CGAL::Object)), 
	    this, SLOT(get_object(CGAL::Object)));
  }
  ~My_window(){delete widget;}
private slots:

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
  My_layer v;                   //an instance of a layer
  CGAL::Qt_widget_standard_toolbar *std_toolbar; 
                                //the standard toolbar
  CGAL::Qt_widget_get_point<Rep> get_point;   
			        //the generic tool that creates Cgal points
  QToolButton *get_point_button;//the toolbar button
};

#include "tutorial8.moc"

int main( int argc, char **argv )
{
    QApplication app( argc, argv );
    My_window W(600,600);
    app.setMainWidget( &W );
    W.show();
    W.setCaption("Using the Standard Toolbar");
    return app.exec();
}
