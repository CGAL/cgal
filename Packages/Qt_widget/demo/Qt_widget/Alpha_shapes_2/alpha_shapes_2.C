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


//STL 
#include <fstream>
#include <stack>
#include <set>
#include <string>

//CGAL
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/squared_distance_2.h> 
#include <CGAL/Point_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/point_generators_2.h>

//Qt_widget
#include "Qt_widget_toolbar_layers.h"
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_widget_Alpha_shape_2.h>
#include "Qt_widget_toolbar.h"



//Qt
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
#include <qslider.h>
#include <qlayout.h>
#include <qinputdialog.h>




typedef double Coord_type;
typedef CGAL::Cartesian<Coord_type>  Rep;
typedef CGAL::Point_2<Rep>  Point;
typedef CGAL::Segment_2<Rep>  Segment;
typedef CGAL::Line_2<Rep>  Line;
typedef CGAL::Triangle_2<Rep>  Triangle;

typedef CGAL::Triangulation_2<Rep> Triangulation;
typedef std::list<Point>     CGALPointlist;



typedef CGAL::Alpha_shape_2<Delaunay>  Alpha_shape;
typedef Alpha_shape::Face  Face;
typedef Alpha_shape::Vertex Vertex;
typedef Alpha_shape::Edge Edge;
typedef Alpha_shape::Face_handle  Face_handle;
typedef Alpha_shape::Vertex_handle Vertex_handle;
typedef Alpha_shape::Face_circulator  Face_circulator;
typedef Alpha_shape::Vertex_circulator  Vertex_circulator;
typedef Alpha_shape::Locate_type Locate_type;
typedef Alpha_shape::Face_iterator  Face_iterator;
typedef Alpha_shape::Vertex_iterator  Vertex_iterator;
typedef Alpha_shape::Edge_iterator  Edge_iterator;
typedef Alpha_shape::Edge_circulator  Edge_circulator;
typedef Alpha_shape::Coord_type Coord_type;
typedef Alpha_shape::Alpha_iterator Alpha_iterator;



const QString my_title_string("Alpha_shapes_2 Demo with"
			      " CGAL Qt_widget");

Delaunay      tr1;
CGALPointlist L;
Alpha_shape   A;
int           current_state;
double        alpha_index;

class Qt_layer_show_alpha_shape : public CGAL::Qt_widget_layer
{
public:
  Qt_layer_show_alpha_shape(){};
private:
  void draw(){
    A.set_alpha(alpha_index);
    A.set_mode(Alpha_shape::GENERAL);
    *widget << CGAL::LineWidth(2) << CGAL::GREEN;
    *widget << A;
  }
};

class MyWindow : public QMainWindow
{
  Q_OBJECT
public:
  MyWindow(int w, int h): win(this) {
  setCentralWidget(&win);

  //create a timer for checking if somthing changed
  QTimer *timer = new QTimer( this );
  connect( timer, SIGNAL(timeout()),
           this, SLOT(timerDone()) );
  timer->start( 200, FALSE );


  // file menu
  QPopupMenu * file = new QPopupMenu( this );
  menuBar()->insertItem( "&File", file );
  file->insertItem("&New", this, SLOT(new_instance()), CTRL+Key_N);
  file->insertItem("New &Window", this, SLOT(new_window()), CTRL+Key_W);
  file->insertSeparator();
  file->insertItem("&Load Triangulation", this, SLOT(load_triangulation()), CTRL+Key_L);
  file->insertItem("&Save Triangulation", this, SLOT(save_triangulation()), CTRL+Key_T);
  file->insertSeparator();
  file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_X );
  file->insertItem( "&Quit", qApp, SLOT( closeAllWindows() ), CTRL+Key_Q );


  // drawing menu
  QPopupMenu * edit = new QPopupMenu( this );
  menuBar()->insertItem( "&Edit", edit );
  edit->insertItem("&Generate Triangulation", this, SLOT(gen_tr()), CTRL+Key_G );
  edit->insertItem("&Change Alpha", this, SLOT(change_alpha()), CTRL+Key_C );

  // help menu
  QPopupMenu * help = new QPopupMenu( this );
  menuBar()->insertItem( "&Help", help );
  help->insertItem("&About", this, SLOT(about()), CTRL+Key_A );
  help->insertItem("About &Qt", this, SLOT(aboutQt()) );

  //the new tools toolbar
  setUsesBigPixmaps(TRUE);
  newtoolbar = new CGAL::Tools_toolbar(&win, this, &tr1);	
  //the new scenes toolbar
  vtoolbar = new CGAL::Layers_toolbar(&win, this, &tr1);
  //the standard toolbar
  stoolbar = new CGAL::Qt_widget_standard_toolbar (&win, this);
  this->addToolBar(stoolbar->toolbar(), Top, FALSE);
  this->addToolBar(newtoolbar->toolbar(), Top, FALSE);
  this->addToolBar(vtoolbar->toolbar(), Top, FALSE);
  
  win << CGAL::LineWidth(2) << CGAL::BackgroundColor (CGAL::BLACK);

  resize(w,h);
  win.show();

  win.setMouseTracking(TRUE);
	
  //connect the widget to the main function that receives the objects
  connect(&win, SIGNAL(new_cgal_object(CGAL::Object)), 
    this, SLOT(get_new_object(CGAL::Object)));

  //application flag stuff
  old_state = 0;
  //layers
  win.attach(&alpha_shape_layer);
  alpha_index = 0.001;
  };

  ~MyWindow()
  {
  };

private:
  void something_changed(){current_state++;};
signals:
  void was_repainted();
  
public slots:
  void set_window(double xmin, double xmax, double ymin, double ymax)
  {
    win.set_window(xmin, xmax, ymin, ymax);
  }
  void new_instance()
  {
    win.lock();
    win.clear();
    tr1.clear();
    A.clear();
    L.clear();
    win.clear_history();
    win.set_window(-1.1, 1.1, -1.1, 1.1); // set the Visible Area to the Interval
    win.unlock();
    something_changed();
  }

private slots:
  void get_new_object(CGAL::Object obj)
  {
    Point p;
    if(CGAL::assign(p,obj)) {
      tr1.insert(p);
      L.push_back(p);
      A.make_alpha_shape(L.begin(), L.end());
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

  void change_alpha() {
    bool ok = FALSE;
    double res = QInputDialog::getDouble( tr( "Please enter a decimal number" ),"Between 0 and 1", 0.001, 0, 1, 3, &ok, this );
    if ( ok ){
      alpha_index = res;
      win.redraw();
    }
  }

  void timerDone()
  {
    if(old_state!=current_state){
      win.redraw();
      old_state = current_state;
    }
  }	

  void gen_tr()
  {
    tr1.clear();
    A.clear();
    win.clear_history();
    win.lock();
    win.set_window(-1.1, 1.1, -1.1, 1.1); // set the Visible Area to the Interval

    // send resizeEvent only on show.
    win.unlock();
    CGAL::Random_points_in_disc_2<Point> g(0.5);
    for(int count=0; count<200; count++) {
      tr1.insert(*g);
      L.push_back(*g++);
    }
    A.make_alpha_shape(L.begin(), L.end());
    win.redraw();
    something_changed();
  }
	
  void save_triangulation()
  {
    QString fileName = 
      QFileDialog::getSaveFileName( "triangulation.cgal", 
				    "Cgal files (*.cgal)", this );
    if ( !fileName.isNull() ) {                 // got a file name
      std::ofstream out(fileName);
      CGAL::set_ascii_mode(out);
      out << tr1 << std::endl;    
    }
  }

	

  void load_triangulation()
  {
    QString s( QFileDialog::getOpenFileName( QString::null,
			    "CGAL files (*.cgal)", this ) );
    if ( s.isEmpty() )
        return;
    tr1.clear();
    std::ifstream in(s);
    CGAL::set_ascii_mode(in);
    in >> tr1;
    something_changed();
  }

	

private:
  
  CGAL::Qt_widget	    win;		
  CGAL::Tools_toolbar	    *newtoolbar;
  CGAL::Layers_toolbar	    *vtoolbar;
  CGAL::Qt_widget_standard_toolbar    *stoolbar;
  //if a CGAL::Point is received should be true
  int			    old_state;
  Qt_layer_show_alpha_shape alpha_shape_layer;		
};

#include "alpha_shapes_2.moc"


int
main(int argc, char **argv)
{
  QApplication app( argc, argv );
    app.setStyle( new QPlatinumStyle );
    QPalette p( QColor( 250, 215, 100 ) );
    app.setPalette( p, TRUE );
  MyWindow win(600, 600); // physical window size
  app.setMainWidget(&win);
  win.setCaption(my_title_string);
  win.setMouseTracking(TRUE);
  win.show();
  // because Qt send resizeEvent only on show.
  win.set_window(-1, 1, -1, 1);
  current_state = -1;
  return app.exec();
}

#endif // CGAL_USE_QT
