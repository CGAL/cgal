// if QT is not installed, a message will be issued in runtime.
#ifndef CGAL_USE_QT
#include <iostream>


int main(int, char*)
{

  std::cout << "Sorry, this demo needs QT...";
  std::cout << std::endl;

  return 0;
}

#else

#include <fstream>
#include <stack>
#include <set>
#include <string>

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/point_generators_2.h>


#include <CGAL/IO/Qt_widget.h>
#include "Qt_widget_toolbar.h"
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_widget_layer.h>

#include <qplatinumstyle.h>
#include <qapplication.h>
#include <qmainwindow.h>
#include <qstatusbar.h>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qfiledialog.h>
#include <qtimer.h>

typedef double Coord_type;
typedef CGAL::Cartesian<Coord_type>  Rep;

typedef CGAL::Point_2<Rep>    Point;
typedef CGAL::Segment_2<Rep>  Segment;

const QString my_title_string("Convex_Hull_2 Demo with"
			      " CGAL Qt_widget");

//global flags and variables
int current_state;
std::list<Point>	  list_of_points;


class Qt_layer_show_ch : public CGAL::Qt_widget_layer
{
public:
  void draw()
  {
     widget->lock();
    *widget << CGAL::PointSize(7) << CGAL::PointStyle(CGAL::CROSS);
    *widget << CGAL::GREEN;
    std::list<Point>::iterator itp = list_of_points.begin();
    while(itp!=list_of_points.end())
    {
      *widget << (*itp++);
    }

    std::list<Point>	out;
    std::list<Segment>	Sl;
    CGAL::convex_hull_points_2(list_of_points.begin(), list_of_points.end(), std::back_inserter(out));

    if( out.size() > 1 ) {
      Point pakt,prev,pstart;

      std::list<Point>::const_iterator it;
      it=out.begin();
      prev= *it; pstart=prev;
      it++;

      for(; it != out.end(); ++it) {
      	pakt= *it;
	      Sl.push_back(Segment(prev,pakt));
	      prev=pakt;
      }
      Sl.push_back(Segment(pakt,pstart));

      *widget << CGAL::RED;
      std::list<Segment>::iterator its = Sl.begin();
      while(its!=Sl.end())
      {
        *widget << (*its++);
      }
    }
    widget->unlock();
  };	
  
};//end class 

class MyWindow : public QMainWindow
{
  Q_OBJECT
public:
  MyWindow(int w, int h){
    widget = new CGAL::Qt_widget(this);
    setCentralWidget(widget);
    
    //create a timer for checking if somthing changed
    QTimer *timer = new QTimer( this );
    connect( timer, SIGNAL(timeout()),
           this, SLOT(timer_done()) );
    timer->start( 200, FALSE );

    // file menu
    QPopupMenu * file = new QPopupMenu( this );
    menuBar()->insertItem( "&File", file );
    file->insertItem("&New", this, SLOT(new_instance()), CTRL+Key_N);
    file->insertItem("New &Window", this, SLOT(new_window()), CTRL+Key_W);
    file->insertSeparator();
    file->insertItem("Print", widget, SLOT(print_to_ps()), CTRL+Key_P);
    file->insertSeparator();
    file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_X );
    file->insertItem( "&Quit", qApp, SLOT( closeAllWindows() ), CTRL+Key_Q );

    // drawing menu
    QPopupMenu * draw = new QPopupMenu( this );
    menuBar()->insertItem( "&Draw", draw );
    draw->insertItem("&Generate points", this, SLOT(gen_points()), CTRL+Key_G );

    // help menu
    QPopupMenu * help = new QPopupMenu( this );
    menuBar()->insertItem( "&Help", help );
    help->insertItem("&About", this, SLOT(about()), CTRL+Key_A );
    help->insertItem("About &Qt", this, SLOT(aboutQt()) );

    //the new tools toolbar
    newtoolbar = new CGAL::Tools_toolbar(widget, this, &list_of_points);	
    //the standard toolbar
    stoolbar = new CGAL::Qt_widget_standard_toolbar (widget, this);
    this->addToolBar(stoolbar->toolbar(), Top, FALSE);
    this->addToolBar(newtoolbar->toolbar(), Top, FALSE);
  
    *widget << CGAL::LineWidth(2) << CGAL::BackgroundColor (CGAL::BLACK);
  
    resize(w,h);

    widget->setMouseTracking(TRUE);
	
    //connect the widget to the main function that receives the objects
    connect(widget, SIGNAL(new_cgal_object(CGAL::Object)), 
      this, SLOT(get_new_object(CGAL::Object)));

    //application flag stuff
    old_state = 0;

    //layers
    widget->attach(&testlayer);
  };

private:
  void something_changed(){current_state++;};
  
public slots:
  void set_window(double xmin, double xmax, double ymin, double ymax)
  {
    widget->set_window(xmin, xmax, ymin, ymax);
  }
  void new_instance()
  {
    widget->lock();
    list_of_points.clear();
    widget->set_window(-1.1, 1.1, -1.1, 1.1); // set the Visible Area to the Interval
    widget->unlock();
    something_changed();
  }

private slots:
  void get_new_object(CGAL::Object obj)
  {
    Point p;
    if(CGAL::assign(p,obj)) {
      list_of_points.push_back(p);
      something_changed();
    }
  };

  void about()
  {
    QMessageBox::about( this, my_title_string,
		"This is a demo for Triangulation,\n"
  		"Copyright CGAL @2001");
  };

  void aboutQt()
  {
    QMessageBox::aboutQt( this, my_title_string );
  }

  void new_window(){
    MyWindow *ed = new MyWindow(500, 500);
    ed->setCaption("Layer");
    ed->show();
    ed->set_window(-1.1, 1.1, -1.1, 1.1);
    something_changed();
  }

  void timer_done()
  {
    if(old_state!=current_state){
      widget->redraw();
      old_state = current_state;
    }
  }	

  void gen_points()
  {
    widget->clear_history();
    widget->set_window(-1.1, 1.1, -1.1, 1.1); // set the Visible Area to the Interval

    // send resizeEvent only on show.
    CGAL::Random_points_in_disc_2<Point> g(0.5);
    for(int count=0; count<200; count++) {
      list_of_points.push_back(*g++);
    }
    something_changed();
  }
	
	

private:
  CGAL::Qt_widget					*widget;
  CGAL::Tools_toolbar				*newtoolbar;
  CGAL::Qt_widget_standard_toolbar  *stoolbar;
  int								old_state;
  Qt_layer_show_ch					testlayer;
};

#include "convex_hull_2.moc"


int
main(int argc, char **argv)
{
  QApplication app( argc, argv );
    app.setStyle( new QPlatinumStyle );
    QPalette p( QColor( 250, 215, 100 ) );
    app.setPalette( p, TRUE );
  MyWindow widget(500,500); // physical window size
  app.setMainWidget(&widget);
  widget.setCaption(my_title_string);
  widget.setMouseTracking(TRUE);
  widget.show();
  // because Qt send resizeEvent only on show.
  widget.set_window(-1, 1, -1, 1);
  current_state = -1;
  return app.exec();
}

#endif // CGAL_USE_QT
