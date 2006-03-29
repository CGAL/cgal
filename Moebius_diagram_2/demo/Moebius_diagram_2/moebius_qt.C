#include <iostream>
#include <fstream>

#include <CGAL/basic.h>

// if QT is not installed, a message will be issued in runtime.
#ifndef CGAL_USE_QT


int main(int, char*)
{

  std::cout << "Sorry, this demo needs QT...";
  std::cout << std::endl;

  return 0;
}

#else

#include "moebius_cgal_types.h"
#include "timer.h"

//Qt_widget headers
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
//#include <CGAL/IO/Qt_widget_helpwindow.h>

//#include "regular_triangulation_2_toolbar.h"
#include "moebius_2_qt_toolbar.h"
//#include "regular_triangulation_2_toolbar_layers.h"
#include "moebius_2_qt_toolbar_layers.h"

//Qt headers
//#include <qplatinumstyle.h>
#include <qcommonstyle.h>

#include <qapplication.h>
#include <qmainwindow.h>
#include <qstatusbar.h>
#include <qmessagebox.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qfiledialog.h>
#include <qtimer.h>


#include "psstream.h"

const QString title_string("Moebius_diagram_2");

MD md;
int		 current_state;
double xmin, ymin, xmax, ymax;


void print_delay (std::ostream &os, double delay, struct tms *_tms)
{
  /*
  os << timer::convert (_tms->tms_utime) << " user, "
     << timer::convert (_tms->tms_stime) << " sys, "
     << delay << " total";
  */
}

bool bounding_box (double &xmin, double &ymin, double &xmax, double ymax)
{
  if (! md.initialized()) return false;

  Vertex_iterator it = md.rt().vertices_begin();

  xmin = xmax = CGAL::to_double ((*it).point().x());
  ymin = ymax = CGAL::to_double ((*it).point().y());
  Vertex_iterator vend = md.rt().vertices_end();
  while(it != vend) {
    if (md.is_finite (it)) {
      if(xmin > CGAL::to_double ((*it).point().x()))
	xmin = CGAL::to_double ((*it).point().x());
      if(xmax < CGAL::to_double ((*it).point().x()))
	xmax = CGAL::to_double ((*it).point().x());
      if(ymin > CGAL::to_double ((*it).point().y()))
	ymin = CGAL::to_double ((*it).point().y());
      if(ymax < CGAL::to_double ((*it).point().y()))
	ymax = CGAL::to_double ((*it).point().y());
    }
    it++;
  }
  return true;
}

void load (std::ifstream &fin)
{
  std::istream_iterator<WPoint> start (fin);
  std::istream_iterator<WPoint> stop;

  struct tms _tms;
  double _delay;
  timer _timer, _total;

  _total.restart ();
  _timer.restart ();
  md.init (start, stop);
  _delay = _timer.elapsed (&_tms);
  std::cout << "     init: ";
  print_delay (std::cout, _delay, &_tms);
  std::cout << "\n";
  
  _timer.restart ();
  md.build ();
  _delay = _timer.elapsed (&_tms);
  std::cout << "    build: ";
  print_delay (std::cout, _delay, &_tms);
  std::cout << "\n";

  _timer.restart ();
  md.construct ();
  _delay = _timer.elapsed (&_tms);
  std::cout << "construct: ";
  print_delay (std::cout, _delay, &_tms);
  std::cout << "\n";

  _delay = _total.elapsed (&_tms);
  std::cout << "    total: ";
  print_delay (std::cout, _delay, &_tms);
  std::cout << "\n";
}

class Window : public QMainWindow
{
  Q_OBJECT
public:
  Window(int w, int h, int lw = 1, int ps = 4, CGAL::PointStyle pt = CGAL::PLUS)
    { 
      widget = new CGAL::Qt_widget(this);
      setCentralWidget(widget);

      widget->setLineWidth (lw);
      widget->setPointSize (ps);
      widget->setPointStyle (pt);

      // file menu
      QPopupMenu * file = new QPopupMenu( this );
      menuBar()->insertItem( "&File", file );
      file->insertItem("&New", this, SLOT(new_instance()), CTRL+Key_N);
      file->insertItem("New &Window", this, SLOT(new_window()), CTRL+Key_W);
      file->insertSeparator();
      file->insertItem("&Load Points", this, 
		       SLOT(load_triangulation()), CTRL+Key_L);
      //    file->insertItem("&Save Triangulation", this,
      //		      SLOT(save_triangulation()), CTRL+Key_S);
      file->insertSeparator();
      file->insertItem("Export to EPS", this, SLOT(my_export()), CTRL+Key_E);
      file->insertItem("Print", widget, SLOT(print_to_ps()), CTRL+Key_P);
      file->insertSeparator();
      file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_X );
      file->insertItem( "&Quit", qApp,
			SLOT( closeAllWindows() ), CTRL+Key_Q );

      // help menu
      QPopupMenu * help = new QPopupMenu( this );
      menuBar()->insertItem( "&Help", help );
      help->insertItem("How To", this, SLOT(howto()), Key_F1);
      help->insertSeparator();
      help->insertItem("&About", this, SLOT(about()), CTRL+Key_A );
      help->insertItem("About &Qt", this, SLOT(aboutQt()) );

      *widget << CGAL::BackgroundColor(CGAL::BLACK);
      resize(w, h);
      widget->set_window(-1, 1, -1, 1);
      widget->setMouseTracking(TRUE);

      //the standard toolbar
      stoolbar = new CGAL::Qt_widget_standard_toolbar (widget, this, "ST");
      //the input layers toolbar
      //newtoolbar = new Tools_toolbar(widget, this, &md);
      //the other layers toolbar
      ltoolbar = new Layers_toolbar(widget, this, &md);


      //connect the widget to the main function that receives the objects
      connect(widget, SIGNAL(new_cgal_object(CGAL::Object)), 
	      this, SLOT(get_new_object(CGAL::Object)));
    
      //create a timer to check if something changed
      QTimer *timer = new QTimer( this );
      connect( timer, SIGNAL(timeout()),
	       this, SLOT(timerDone()) );
      timer->start( 200, FALSE );
      old_state = 0; 

    }
  void 
  init_coordinates(){
    xmin = -1; xmax = 1;
    ymin = -1; ymax = 1;
    bounding_box (xmin, ymin, xmax, ymax);
  }
  private slots:
  void 
  new_instance(){
    widget->lock();
    widget->clear();
    widget->clear_history();
    md.clear();
    // set the Visible Area to the Interval
    widget->set_window(-1.1, 1.1, -1.1, 1.1);
    widget->unlock();
    something_changed();
  }
  void 
  new_window(){
    Window *ed = new Window(500, 500);
    ed->setCaption("Layer");
    bounding_box (xmin, ymin, xmax, ymax);
    ed->widget->clear_history();
    ed->widget->set_window(xmin, xmax, ymin, ymax);
    ed->show();
    something_changed();
  }
  
  void 
  timerDone(){
    if(old_state!=current_state){
      widget->redraw();
      old_state = current_state;
    }
  }

  void 
  howto(){
    //     QString home;
    //     home = "help/rindex.html";
    //     HelpWindow *help = new HelpWindow(home, ".", 0, "help viewer");
    //     help->resize(400, 400);
    //     help->setCaption("Demo HowTo");
    //     help->show();
  }

  void 
  about(){
    QMessageBox::about( this, title_string,
			"This is a demo for Moebius_diagram_2\n"
			"Copyright NOBODY @2003");
}

void 
aboutQt(){
  QMessageBox::aboutQt( this, title_string );
}

void 
get_new_object(CGAL::Object obj){
#if 0
  Point p;
  Circle c;
  if(CGAL::assign(p,obj)) {
    md.insert(p);
    something_changed();
  } else if (CGAL::assign(c, obj)){
    md.insert(Gt::Weighted_point(Point(c.center()), c.squared_radius()));
    something_changed();
  }
#endif
}

  
  void
  my_export(){
    std::cout << "et paf !\n";
    QString s( QFileDialog::getOpenFileName( QString::null,
					     "EPS files (*.eps)", this ) );
    if ( s.isEmpty() )
      return;
    const char *file = s.ascii ();
    if (! md.initialized()) return;
    psstream ps (file, xmin, ymin, xmax, ymax);
    {
      Vertex_iterator i = md.rt().vertices_begin (), end = md.rt().vertices_end ();
      while (i != end) {
	if (md.is_finite (i))
	  ps << (Point) (*i).point();
	++i;
      }
    }
    if (! md.constructed()) return;
    {
      //      Edge_iterator i = md.edges_begin (), end = md.edges_end ();
      //      while (i != end) {
	
      // }
    }
  }

  void
    load_triangulation(){
    QString s( QFileDialog::getOpenFileName( QString::null,
					     "CGAL files (*.cgal)", this ) );
    if ( s.isEmpty() )
      return;
    md.clear();

    const char *file = s.ascii ();
    std::ifstream fin (file);
    std::cout << "loading " << file << "...\n";
    load (fin);
    bounding_box (xmin, ymin, xmax, ymax);
    something_changed ();
  }
  void
    save_triangulation(){
  }
 private:
  inline  void something_changed(){current_state++;};
  CGAL::Qt_widget *widget;
  CGAL::Qt_widget_standard_toolbar  *stoolbar;
  Tools_toolbar                     *newtoolbar;
  Layers_toolbar                    *ltoolbar;
  int                               old_state;

};

#include "moebius_qt.moc"

int
main(int argc, char **argv)
{
  int linewidth = 1, pointsize = 4;
  CGAL::PointStyle style = CGAL::PLUS;

  if (argc > 2 && !strcmp (argv[1], "-style")) {
    std::cout << "parsing style string `"<<argv[2]<<"'";
    char *pointstyle;
    sscanf (strtok (argv[2], ":"), "%d", &linewidth);
    std::cout << "  line width = " << linewidth << "\n";
    sscanf (strtok (NULL, ":"), "%d", &pointsize);
    std::cout << "  point size = " << pointsize << "\n";
    pointstyle = strtok (NULL, ":");
    if (!strcmp (pointstyle, "cross")) style = CGAL::CROSS;
    else if (!strcmp (pointstyle, "plus")) style = CGAL::PLUS;
    std::cout << "  point style = " << pointstyle << " (" << style << ")\n" ;

    argc -= 2; argv += 2;
  }
 if (argc > 1) {
    std::cout << "loading " << argv[1] << "...\n";
    std::ifstream fin (argv[1]);
    load (fin);
    argc--; argv++;
  }

  std::cout << "Initializing Qt window...\n";
  QApplication app( argc, argv );

  //  std::cout << "Window W()...\n";
  Window W(500,500); // physical widget size

  //std::cout << "app.setMainWidget()...\n";
  app.setMainWidget(&W);

  
  //std::cout << "W.setCaption(title_string)...\n";
  W.setCaption(title_string);

  //std::cout << "W.setMouseTracking(TRUE)...\n";
  W.setMouseTracking(TRUE);

  //std::cout << "W.show(); ...\n";
  W.show();  

  //std::cout << "W.init_coordinates(); ...\n";
  W.init_coordinates();  

  //std::cout << "current_state = -1...\n";
  current_state = -1;

  std::cout << "done.\n";
  return app.exec();
}


#endif //CGAL_USE_QT
