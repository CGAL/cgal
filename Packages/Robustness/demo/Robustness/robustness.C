// file          : robustness.C
// author(s)     : Stefan Schirra, Radu Ursu

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
#include "cgal_types.h"
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_help_window.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/pixmaps/demoicon.xpm>

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
#include <qprogressdialog.h>


#include <fstream>
#include <stack>
#include <set>
#include <string>

const QString my_title_string("Robustness Demo with"
			      " CGAL Qt_widget");

//global flags and variables
int                           current_state;
std::vector<double_Segment>   double_segments, double_segments2;
std::vector<real_Segment>     real_segments;

std::vector<CartesianFloat::Segment_2>    float_segments, float_segments2;
std::vector<HomogeneousFloat::Segment_2>  hfloat_segments, hfloat_segments2;
std::vector<HomogeneousDouble::Segment_2> hdouble_segments, hdouble_segments2;

std::vector<double_Point >    double_intersection_points;
std::vector<real_Point >      real_intersection_points;
std::vector<double_Point >    double_convex_hull;
std::vector<real_Point >      real_convex_hull;

class show_intersection_points : public CGAL::Qt_widget_layer{
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
    *widget << CGAL::LineWidth(1) << CGAL::BLUE;
    std::vector<double_Segment>::iterator dit2 =
      double_segments2.begin();
    while(dit2!=double_segments2.end()){
      *widget << (*dit2);
      dit2++;
    }
    
    widget->unlock();
  }
};

class show_hull_of_intersection_points : public CGAL::Qt_widget_layer{
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
    draw->insertItem("Convex &Hull of intersection points", this, 
		     SLOT(hull_points()), CTRL+Key_H );
    
    draw->insertItem("&Intersection points on segments", this, 
		     SLOT(intersection_points()), CTRL+Key_I );

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
    widget->attach(&hull_layer);
    widget->attach(&int_points_layer);

    qte = new QTextBrowser(NULL, "INFO");
    qte->setCaption("Information Window");

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
    CGAL::Qt_help_window *help = new 
      CGAL::Qt_help_window(home, ".", 0, "help viewer");
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

  void intersection_points(){
    stoolbar->clear_history();
    widget->set_window(-1.1, 1.1, -1.1, 1.1);
    gen_segments();

    //computing using Cartesian<float>
    CartesianFloat::Point_2      p;
    CartesianFloat::Intersect_2  intersection = 
      CartesianFloat().intersect_2_object();
    CGAL::Timer watch;
    int is_count = 0;
    int ol_count = 0;
    int bl_count = 0;
    watch.start();    
    for( std::vector<CartesianFloat::Segment_2>::const_iterator i 
	   = float_segments.begin(); i != float_segments.end(); ++i)
        for( std::vector<CartesianFloat::Segment_2>::const_iterator j 
	       = float_segments2.begin(); j != float_segments2.end(); ++j)
        {
            if ( CGAL::assign(p, intersection(*i,*j)) )
            {
                is_count++;
                int ok1 = ( i->has_on(p)) ? 1 : 0;
                int ok2 = ( j->has_on(p)) ? 1 : 0;
                bl_count += ok1*ok2;
                ol_count += ok1+ok2;
            }
        }
    watch.stop();
    QString s("%1"), s1("%1");
    qte->resize(400, 300);
    qte->show();
    qte->setText("INFORMATION Cartesian<float>:");
    s.setNum(is_count);
    qte->append("Intersection points found: " + s);
    s.setNum(bl_count); s1.setNum((double)bl_count/is_count * 100);
    qte->append(s + " of them (" + s1 + "%) lie on both segments.");
    s.setNum(2*is_count); s1.setNum(ol_count);
    qte->append("Out of the " + s + " point-on-segment tests, ");
    s.setNum((double)ol_count/is_count * 50);
    qte->append(s1 + "(" + s + "%) are positive.");
    s.setNum(watch.time());
    qte->append("Computation time = " + s + " sec.");


    //computing using Cartesian<double>
    CartesianDouble::Intersect_2  double_intersection = 
      CartesianDouble().intersect_2_object();
    double_Point pd;
    is_count = 0;
    ol_count = 0;
    bl_count = 0;
    watch.reset();
    watch.start();    
    for( std::vector<CartesianDouble::Segment_2>::const_iterator i 
	   = double_segments.begin(); i != double_segments.end(); ++i)
        for( std::vector<CartesianDouble::Segment_2>::const_iterator j 
	       = double_segments2.begin(); j != double_segments2.end(); ++j)
        {
            if ( CGAL::assign(pd, double_intersection(*i,*j)) )
            {
                is_count++;
                int ok1 = ( i->has_on(pd)) ? 1 : 0;
                int ok2 = ( j->has_on(pd)) ? 1 : 0;
                bl_count += ok1*ok2;
                ol_count += ok1+ok2;
            }
        }
    watch.stop();
    qte->append("_____");
    qte->append("INFORMATION Cartesian " + QString::null + "double>:");
    s.setNum(is_count);
    qte->append("Intersection points found: " + s);
    s.setNum(bl_count); s1.setNum((double)bl_count/is_count * 100);
    qte->append(s + " of them (" + s1 + "%) lie on both segments.");
    s.setNum(2*is_count); s1.setNum(ol_count);
    qte->append("Out of the " + s + " point-on-segment tests, ");
    s.setNum((double)ol_count/is_count * 50);
    qte->append(s1 + "(" + s + "%) are positive.");
    s.setNum(watch.time());
    qte->append("Computation time = " + s + " sec.");


    //computing using Homogeneous<float>
    HomogeneousFloat::Intersect_2  hfloat_intersection = 
      HomogeneousFloat().intersect_2_object();
    HomogeneousFloat::Point_2 phf;
    is_count = 0;
    ol_count = 0;
    bl_count = 0;
    watch.reset();
    watch.start();    
    for( std::vector<HomogeneousFloat::Segment_2>::const_iterator i 
	   = hfloat_segments.begin(); i != hfloat_segments.end(); ++i)
        for( std::vector<HomogeneousFloat::Segment_2>::const_iterator j 
	       = hfloat_segments2.begin(); j != hfloat_segments2.end(); ++j)
        {
            if ( CGAL::assign(phf, hfloat_intersection(*i,*j)) )
            {
                is_count++;
                int ok1 = ( i->has_on(phf)) ? 1 : 0;
                int ok2 = ( j->has_on(phf)) ? 1 : 0;
                bl_count += ok1*ok2;
                ol_count += ok1+ok2;
            }
	}
    watch.stop();
    qte->append("_____");
    qte->append("INFORMATION Homogeneous float>:");
    s.setNum(is_count);
    qte->append("Intersection points found: " + s);
    s.setNum(bl_count); s1.setNum((double)bl_count/is_count * 100);
    qte->append(s + " of them (" + s1 + "%) lie on both segments.");
    s.setNum(2*is_count); s1.setNum(ol_count);
    qte->append("Out of the " + s + " point-on-segment tests, ");
    s.setNum((double)ol_count/is_count * 50);
    qte->append(s1 + "(" + s + "%) are positive.");
    s.setNum(watch.time());
    qte->append("Computation time = " + s + " sec.");


    hull_layer.deactivate();
    int_points_layer.activate();
    something_changed();
  }

  void hull_points(){
    stoolbar->clear_history();
    widget->set_window(-1.1, 1.1, -1.1, 1.1);
    QString s("%1");
    qte->resize(400, 300);
    qte->show();
    qte->setText("INFORMATION:");
    qte->append("We compute the intersection points of segments using exact"
                " and double arithmetic. Then we compute the convex hull of"
                " those. If there are points that are different in different"
                " arithmetic, we mark them by a red circle.");
    
    QProgressDialog progress( "Generating segments...", 
      "Cancel computing", 60, NULL, "Compute random segments ...", true );
    progress.setCaption("Progress bar");
    progress.resize(200, 50);
    progress.setTotalSteps(60);
    progress.setProgress(10);
    progress.setMinimumDuration(0);
    progress.setLabelText("Compute random segments...");

    gen_segments();
    progress.setLabelText("Compute intersection points...");
    progress.setProgress(20);
    double_convex_hull.clear();
    real_convex_hull.clear();
    double_intersection_points.clear();
    real_intersection_points.clear();
    CGAL::segment_intersection_points_2(
          double_segments.begin(),
          double_segments.end(),
          std::back_inserter( double_intersection_points),
          C_double() );
    progress.setLabelText("Compute intersection points exact...");
    progress.setProgress(30);
    CGAL::segment_intersection_points_2(
          real_segments.begin(),
          real_segments.end(),
          std::back_inserter( real_intersection_points),
          C_real() );
    progress.setLabelText("Compute hull of intersection points...");
    progress.setProgress(40);
    CGAL::convex_hull_points_2(
          double_intersection_points.begin(),
          double_intersection_points.end(),
          std::back_inserter( double_convex_hull));
    progress.setProgress(50);
    progress.setLabelText("Compute hull of intersection points exact...");
    CGAL::convex_hull_points_2(
          real_intersection_points.begin(),
          real_intersection_points.end(),
          std::back_inserter( real_convex_hull));
    progress.setProgress(60);

    hull_layer.activate();
    int_points_layer.deactivate();
    something_changed();
  }

  void gen_segments()
  {
    Source RS(1);
    Segment_iterator g( RS, RS);
    double_segments.clear();
    real_segments.clear();
    double_segments2.clear();
    float_segments.clear();
    float_segments2.clear();
    hfloat_segments.clear();
    hfloat_segments2.clear();
    hdouble_segments.clear();
    hdouble_segments2.clear();

    CGAL::copy_n( g, 100, std::back_inserter( double_segments) );
    CGAL::copy_n( g, 100, std::back_inserter( double_segments2) );

    CGAL::Cartesian_converter<C_double, C_real> converter;
    std::transform( double_segments.begin(),
                    double_segments.end(),
                    std::back_inserter( real_segments),
                    converter );

    CGAL::Cartesian_converter<C_double, CartesianFloat> fconverter;
    std::transform(double_segments.begin(), double_segments.end(),
		    std::back_inserter(float_segments), fconverter);
    std::transform(double_segments2.begin(), double_segments2.end(),
		    std::back_inserter(float_segments2), fconverter);

    CGAL::Cartesian_double_to_Homogeneous< HomogeneousFloat::RT > 
                                                        hfconverter;
    std::transform(double_segments.begin(), double_segments.end(),
		   std::back_inserter(hfloat_segments), hfconverter);
    std::transform(double_segments2.begin(), double_segments2.end(),
		   std::back_inserter(hfloat_segments2), hfconverter);

    CGAL::Cartesian_double_to_Homogeneous<HomogeneousDouble::RT>
                                                        hdconverter;
    std::transform(double_segments.begin(), double_segments.end(),
		   std::back_inserter(hdouble_segments), hdconverter);
    std::transform(double_segments2.begin(), double_segments2.end(),
		   std::back_inserter(hdouble_segments2), hdconverter);


  }	

private:
  CGAL::Qt_widget           *widget;
  CGAL::Qt_widget_standard_toolbar
                            *stoolbar;
  int                       old_state;
  show_hull_of_intersection_points
                            hull_layer;
  show_intersection_points  int_points_layer;
  QTextBrowser              *qte;
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
  QPixmap cgal_icon = QPixmap((const char**)demoicon_xpm);
  widget.setIcon(cgal_icon);
  widget.show();
  current_state = -1;
  return app.exec();
}

#endif // CGAL_USE_QT
