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
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/point_generators_2.h>



#include <CGAL/IO/Qt_Widget.h>
#include "Qt_widget_toolbar.h"
#include "Qt_widget_toolbar_views.h"
#include <CGAL/IO/Qt_widget_standard_toolbar.h>

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

typedef CGAL::Point_2<Rep>  Point;
typedef CGAL::Segment_2<Rep>  Segment;
typedef CGAL::Line_2<Rep>  Line;
typedef CGAL::Triangle_2<Rep>  Triangle;
typedef CGAL::Circle_2<Rep> Circle;

typedef CGAL::Triangulation_2<Rep> Triangulation;
typedef CGAL::Delaunay_triangulation_2<Rep>  Delaunay;



typedef Delaunay::Face_handle               Face_handle;
typedef Delaunay::Vertex_handle             Vertex_handle;
typedef Delaunay::Edge                      Edge;
typedef Triangulation::Line_face_circulator  Line_face_circulator;


const QString my_title_string("Triangulation Demo with"
			      " CGAL Qt_scenes_widget");


Delaunay	tr1;
int		current_state;

class MyWindow : public QMainWindow
{
  Q_OBJECT
public:
  MyWindow(int w, int h): win(this) {
  setCentralWidget(&win);
    
  connect(&win, SIGNAL(mouseReleased(QMouseEvent*)), this,
          SLOT(insert_after_show_conflicts(QMouseEvent*)));

  connect(&win, SIGNAL(resized()), this, SLOT(redrawWin()));
	
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
  QPopupMenu * draw = new QPopupMenu( this );
  menuBar()->insertItem( "&Draw", draw );
  draw->insertItem("&Generate_triangulation", this, SLOT(gen_tr()), CTRL+Key_G );

  // help menu
  QPopupMenu * help = new QPopupMenu( this );
  menuBar()->insertItem( "&Help", help );
  help->insertItem("&About", this, SLOT(about()), CTRL+Key_A );
  help->insertItem("About &Qt", this, SLOT(aboutQt()) );

  //the new tools toolbar
  setUsesBigPixmaps(TRUE);
  newtoolbar = new CGAL::Tools_toolbar(&win, this, &tr1);	
  //the new scenes toolbar
  vtoolbar = new CGAL::Views_toolbar(&win, this, &tr1);
  //the standard toolbar
  stoolbar = new CGAL::Standard_toolbar (&win, this);
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
  got_point = FALSE;
  old_state = 0;
  };

  ~MyWindow()
  {
  };

  void show_conflicts(Point p)
  {
    if(tr1.dimension()<2) return;
    std::list<Face_handle> conflict_faces;
    std::list<Edge>  hole_bd;
    tr1.get_conflicts_and_boundary(p, 
    std::back_inserter(conflict_faces),
    std::back_inserter(hole_bd));
    std::list<Face_handle>::iterator fit = conflict_faces.begin();
    std::list<Edge>::iterator eit = hole_bd.begin();
    win << CGAL::WHITE ;
    for( ; fit != conflict_faces.end(); fit++)  {
      if(! tr1.is_infinite( *fit))
	win << tr1.triangle( *fit );
    }
    win << CGAL::YELLOW;
    for( ; eit != hole_bd.end(); eit++)  {
      if(! tr1.is_infinite( *eit ))
	win << tr1.segment( *eit );
    }
		
  }
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
    win.detach_current_tool();
    win.lock();
    win.clear();
    tr1.clear();
    win.set_window(-1.1, 1.1, -1.1, 1.1); // set the Visible Area to the Interval
    win.unlock();
    something_changed();
  }

  void redrawWin()
  {
    //emit(was_repainted());	
  }
	

private slots:
  void get_new_object(CGAL::Object obj)
  {
    Point p;
    Segment s;
    Line l;
    if (CGAL::assign(l,obj))
    {
      if (tr1.dimension()<1) return;
      win.redraw();
      win.lock();
      Line_face_circulator lfc = 
	tr1.line_walk(l.point(1), l.point(2)), done(lfc);
      if(lfc == (CGAL_NULL_TYPE) NULL){
      } else {
	win << CGAL::BLUE;
	win << CGAL::FillColor(CGAL::WHITE);
	do{
	  if(! tr1.is_infinite( lfc  )){
	    win << tr1.triangle( lfc );
	  }
	}while(++lfc != done);
      }
      win << CGAL::GREEN << l ;
      win << CGAL::noFill;
      win.unlock();
    } else if(CGAL::assign(p,obj)) {
      got_point = TRUE;
      show_conflicts(p);
      tr1.insert(p);
    } 
  };

  void insert_after_show_conflicts(QMouseEvent*)
  {
    if(got_point)
    {
      got_point = FALSE;
      win.redraw();
      something_changed();
    }
  }

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
    ed->setCaption("View");
    ed->show();
    something_changed();
  }

  void timerDone()
  {
    if(old_state!=current_state){
      redrawWin();
      win.redraw();
      old_state = current_state;
    }
  }	

  void gen_tr()
  {
    tr1.clear();
    win.lock();
    win.set_window(-1.1, 1.1, -1.1, 1.1); // set the Visible Area to the Interval

    // send resizeEvent only on show.
    win.unlock();
    CGAL::Random_points_in_disc_2<Point> g(0.5);
    for(int count=0; count<200; count++) {
      tr1.insert(*g++);
    }
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
  
  CGAL::Qt_widget	  win;		
  CGAL::Tools_toolbar	  *newtoolbar;
  CGAL::Views_toolbar	  *vtoolbar;
  CGAL::Standard_toolbar  *stoolbar;
  bool			  got_point;	
	  //if a CGAL::Point is received should be true
  int			  old_state;
	
	
};

#include "triangulationdemo.moc"


int
main(int argc, char **argv)
{
  QApplication app( argc, argv );
    app.setStyle( new QPlatinumStyle );
    QPalette p( QColor( 250, 215, 100 ) );
    app.setPalette( p, TRUE );
  MyWindow win(800,800); // physical window size
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
