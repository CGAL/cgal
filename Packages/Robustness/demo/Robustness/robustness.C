// ============================================================================
//
// Copyright (c) 1997-2003 The CGAL Consortium
// This software and related documentation are part of the Computational
// Geometry Algorithms Library (CGAL).
// This software and documentation are provided "as-is" and without warranty
// of any kind. In no event shall the CGAL Consortium be liable for any
// damage of any kind. 
// ----------------------------------------------------------------------
//
// file          : robustness.C
// package       : Qt_widget
// author(s)     : Radu Ursu
// release_date  : 2002, May 16
//
// coordinator   : Laurent Rineau
//
// email         : contact@cgal.org
// www           : http://www.cgal.org
//
// ======================================================================

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
#  include <CGAL/MP_Float.h>
#  include <CGAL/Quotient.h>
typedef CGAL::Quotient<CGAL::MP_Float> exact_NT;
#include <CGAL/segment_intersection_points_2.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/copy_n.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Cartesian_converter.h> 



#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_widget_helpwindow.h>
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

  typedef CGAL::Cartesian<double>       C_double;
  typedef C_double::Point_2             double_Point;
  typedef C_double::Segment_2           double_Segment;
  typedef CGAL::Cartesian<exact_NT>     C_real;
  typedef C_real::Point_2               real_Point;
  typedef C_real::Segment_2             real_Segment;
  typedef CGAL::Creator_uniform_2<double, double_Point>
                                        Point_creator;
  typedef CGAL::Random_points_in_square_2<double_Point, Point_creator>
                                        Source;
  typedef CGAL::Creator_uniform_2<double_Point,  double_Segment>
                                        Segment_creator;
  typedef CGAL::Join_input_iterator_2<Source, Source, Segment_creator>
                                        Segment_iterator;

const QString my_title_string("Robustness Demo with"
			      " CGAL Qt_widget");

//global flags and variables
int                           current_state;
std::vector<double_Segment>   double_segments;
std::vector<real_Segment>     real_segments;
std::vector<double_Point >    double_intersection_points;
std::vector<real_Point >      real_intersection_points;
std::vector<double_Point >    double_convex_hull;
std::vector<real_Point >      real_convex_hull;

class show_segments : public CGAL::Qt_widget_layer{
public:
  void draw(){
    widget->lock();
    *widget << CGAL::LineWidth(1) << CGAL::GREEN;
    std::vector<double_Segment>::iterator dit =
      double_segments.begin();
    while(dit!=double_segments.end()){
      *widget << (*dit);
      dit++;
    }

    std::list<double_Segment>	Sl;
    if( double_convex_hull.size() > 1 ) {
      double_Point pakt,prev,pstart;

      std::vector<double_Point>::iterator it;
      it=double_convex_hull.begin();
      prev= *it; pstart=prev;
      it++;

      for(; it != double_convex_hull.end(); ++it) {
      	pakt = *it;
        Sl.push_back(double_Segment(prev,pakt));
        prev=pakt;
      }
      Sl.push_back(double_Segment(pakt,pstart));

      *widget << CGAL::BLUE;
      std::list<double_Segment>::iterator its = Sl.begin();
      while(its!=Sl.end())
      {
        *widget << (*its++);
      }
    }

    std::vector<double_Point>::iterator dip =
      double_convex_hull.begin();
    *widget << CGAL::PointStyle(CGAL::CROSS) << CGAL::LineWidth(2);
    *widget << CGAL::WHITE << CGAL::PointSize(8);
    while(dip!=double_convex_hull.end()){
      *widget << (*dip);
      dip++;
    }

    if ( real_convex_hull.size() != double_convex_hull.size() )
    {
      *widget <<CGAL::PointSize(11) << CGAL::PointStyle(CGAL::CIRCLE) 
	      << CGAL::RED;
      std::vector<double_Point >::iterator  dble_it;
      std::vector<real_Point >::iterator    real_it;
      real_it = real_convex_hull.begin();
      for ( dble_it  = double_convex_hull.begin();
          dble_it != double_convex_hull.end();
          ++dble_it )
      {
        if (   (real_it == real_convex_hull.end())
            || (( CGAL::squared_distance(
                      *dble_it,
                      double_Point( CGAL::to_double(real_it->x()),
                                    CGAL::to_double(real_it->y()) ))
               ) > 0.125 )
           )
        {
          *widget << *dble_it;
        } else {
          if ( real_it != real_convex_hull.end() ){
	    ++real_it;
	  }
        }
      } 
    }
    widget->unlock();
  }
};


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
    draw->insertItem("&Generate segments", this, 
		     SLOT(gen_segments()), CTRL+Key_S );

    // help menu
    QPopupMenu * help = new QPopupMenu( this );
    menuBar()->insertItem( "&Help", help );
    help->insertItem("How To", this, SLOT(howto()), Key_F1);
    help->insertSeparator();
    help->insertItem("&About", this, SLOT(about()), CTRL+Key_A );
    help->insertItem("About &Qt", this, SLOT(aboutQt()) );

    //the standard toolbar
    stoolbar = new CGAL::Qt_widget_standard_toolbar (widget, this, "ST");
  
    *widget << CGAL::LineWidth(2) << CGAL::BackgroundColor (CGAL::BLACK);
  
    resize(w,h);
    widget->set_window(-1, 1, -1, 1);
    widget->setMouseTracking(TRUE);
	
    //application flag stuff
    old_state = 0;

    //layers
    widget->attach(&segments_layer);
  };

private:
  void something_changed(){current_state++;};
  
public slots:
  void new_instance()
  {
    widget->lock();
    stoolbar->clear_history();
    double_segments.clear();
    real_segments.clear();
    double_convex_hull.clear();
    real_convex_hull.clear();
    double_intersection_points.clear();
    real_intersection_points.clear();
    widget->set_window(-1.1, 1.1, -1.1, 1.1); 
    // set the Visible Area to the Interval
    widget->unlock();
    something_changed();
  }

private slots:

  void howto(){
    QString home;
    home = "help/index.html";
    HelpWindow *help = new HelpWindow(home, ".", 0, "help viewer");
    help->resize(400, 400);
    help->setCaption("Demo HowTo");
    help->show();
  }

  void about()
  {
    QMessageBox::about( this, my_title_string,
		"This is a demo for Robustness\n"
  		"Copyright CGAL @2003");
  };

  void aboutQt()
  {
    QMessageBox::aboutQt( this, my_title_string );
  }

  void new_window(){
    MyWindow *ed = new MyWindow(500, 500);
    ed->setCaption("Layer");
    ed->stoolbar->clear_history();
    ed->widget->set_window(-1.1, 1.1, -1.1, 1.1);
    ed->show();
    something_changed();
  }

  void timer_done()
  {
    if(old_state!=current_state){
      widget->redraw();
      old_state = current_state;
    }
  }	

  void gen_segments()
  {
    stoolbar->clear_history();
    widget->set_window(-1.1, 1.1, -1.1, 1.1); 
    // set the Visible Area to the Interval
    Source RS(1);
    Segment_iterator g( RS, RS);
    double_segments.clear();
    real_segments.clear();
    double_convex_hull.clear();
    real_convex_hull.clear();
    double_intersection_points.clear();
    real_intersection_points.clear();
    CGAL::copy_n( g, 100, std::back_inserter( double_segments) );
    CGAL::Cartesian_converter<C_double, C_real> converter;
    std::transform( double_segments.begin(),
                    double_segments.end(),
                    std::back_inserter( real_segments),
                    converter );
    CGAL::segment_intersection_points_2(
          double_segments.begin(),
          double_segments.end(),
          std::back_inserter( double_intersection_points),
          C_double() );
    CGAL::segment_intersection_points_2(
          real_segments.begin(),
          real_segments.end(),
          std::back_inserter( real_intersection_points),
          C_real() );
    CGAL::convex_hull_points_2(
          double_intersection_points.begin(),
          double_intersection_points.end(),
          std::back_inserter( double_convex_hull));
    CGAL::convex_hull_points_2(
          real_intersection_points.begin(),
          real_intersection_points.end(),
          std::back_inserter( real_convex_hull));

    something_changed();
  }	

private:
  CGAL::Qt_widget           *widget;
  CGAL::Qt_widget_standard_toolbar
                            *stoolbar;
  int                       old_state;
  show_segments             segments_layer;
};

#include "robustness.moc"


int
main(int argc, char **argv)
{
  QApplication app( argc, argv );
  MyWindow widget(500,500); // physical window size
  app.setMainWidget(&widget);
  widget.setCaption(my_title_string);
  widget.setMouseTracking(TRUE);
  widget.show();
  current_state = -1;
  return app.exec();
}

#endif // CGAL_USE_QT
