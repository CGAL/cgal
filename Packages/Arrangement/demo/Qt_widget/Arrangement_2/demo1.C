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
// file          : main.C
// package       : Qt_widget
// author(s)     : Radu Ursu
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

#include "cgal_types1.h"

#include <fstream>
#include <stack>
#include <set>
#include <string>
#include <list>


#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_Polygon_2.h>
#include "Qt_widget_toolbar1.h"
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_help_window.h>
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



const QString my_title_string("Arrangement Demo with"
			      " CGAL Qt_widget");

//global flags and variables
int                 current_state;
std::list<Curve>    list_of_segments;
Arr                 arr;
bool                pl_valid=false;
Point               pl_point;

class Qt_layer_show_ch : public CGAL::Qt_widget_layer
{
public:
	
  Qt_layer_show_ch(){};


  void draw()
  {
    widget->lock();

    *widget << CGAL::GREEN;
    *widget << CGAL::LineWidth(1);
    std::list<Curve>::iterator itp = list_of_segments.begin();
    while(itp!=list_of_segments.end()){
      *widget << (*itp++);
    }

    if(pl_valid && !(arr.halfedges_begin() == arr.halfedges_end()) ) 
    {
      *widget << CGAL::LineWidth(3);
      *widget << CGAL::YELLOW;

      Arr::Locate_type lt;
      Arr::Point temp_p(pl_point.x(),pl_point.y());
      Arr::Halfedge_handle e = arr.locate(temp_p, lt);
//      std::cout << "locate type " << lt << std::endl;
	
      //color the face on the screen
      Arr::Face_handle f=e->face();

	
      if (f->does_outer_ccb_exist()) 
      {
	Arr::Ccb_halfedge_circulator cc=f->outer_ccb();
	do {
	  *widget << cc->curve();
	} while (++cc != f->outer_ccb());
	
      }

      Arr::Holes_iterator hit=f->holes_begin(),eit=f->holes_end();
      for (;hit!=eit; ++hit) 
      {
	Arr::Ccb_halfedge_circulator cc=*hit; 
	do 
	{
	  *widget << cc->curve();
	} while (++cc != *hit);
      }
      *widget << CGAL::LineWidth(1);
    }

    widget->unlock();
  };	

  void mousePressEvent(QMouseEvent *e)
  {
    if(e->button() == Qt::RightButton)
    {
      if( list_of_segments.empty() )
	return;

      NT x=static_cast<NT>(widget->x_real(e->x()));
      NT y=static_cast<NT>(widget->y_real(e->y()));

      Point p(x,y);
      NT min_dist=100000000;
      std::list<Curve>::iterator itp = list_of_segments.begin();
      std::list<Curve>::iterator it_seg=NULL;
      while(itp!=list_of_segments.end())
      {
	NT dist = CGAL::squared_distance( p, (*itp));
	if( dist < min_dist)
	{
	  min_dist = dist;
	  it_seg = itp;
	}
	itp++;
      }

      
      Arr::Curve_iterator ci = arr.curve_node_begin();
      while(ci != arr.curve_node_end() )
      {
	if( (*ci).curve() == 
	  (*it_seg) )
	{
	  arr.remove_curve( ci );
	  break;
	}
	ci++;
      }

      list_of_segments.erase( it_seg );
      
      pl_valid = false;
      (*widget).redraw();
    }
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
    draw->insertItem("&Generate segments", this,
				SLOT(gen_segments()), CTRL+Key_G );

    // help menu
    QPopupMenu * help = new QPopupMenu( this );
    menuBar()->insertItem( "&Help", help );
    help->insertItem("How To", this, SLOT(howto()), Key_F1);
    help->insertSeparator();
    help->insertItem("&About", this, SLOT(about()), CTRL+Key_A );
    help->insertItem("About &Qt", this, SLOT(aboutQt()) );

    //the standard toolbar
    stoolbar = new CGAL::Qt_widget_standard_toolbar (widget, this, "ST");
    //the new tools toolbar
    newtoolbar = new Tools_toolbar(widget, this, &list_of_segments);	
  
    *widget << CGAL::LineWidth(2) << CGAL::BackgroundColor (CGAL::BLACK);
  
    resize(w,h);
    widget->set_window(-1, 1, -1, 1);
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
  void new_instance()
  {
    widget->lock();
    list_of_segments.clear();
    stoolbar->clear_history();
    widget->set_window(-1.1, 1.1, -1.1, 1.1); 
			// set the Visible Area to the Interval
    widget->unlock();
    something_changed();
  }

private slots:
  void get_new_object(CGAL::Object obj)
  {
    pl_valid = false;

    Segment s;
    if(CGAL::assign(s,obj)) {
      list_of_segments.push_back(s);
      arr.insert(s);
      something_changed();
    }

    Point p;
    if(CGAL::assign(p,obj)) {

      pl_point = p;
      pl_valid = true;
      
      something_changed();
    }

  };

  void about()
  {
    QMessageBox::about( this, my_title_string,
		"This is a demo for the Arrangement package\n"
  		"Copyright CGAL @2003");
  };

  void aboutQt()
  {
    QMessageBox::aboutQt( this, my_title_string );
  }

  void howto(){
    QString home;
    home = "help/index.html";
    CGAL::Qt_help_window * help =
      new CGAL::Qt_help_window(home, ".", 0, "help viewer");
    help->resize(400, 400);
    help->setCaption("Demo HowTo");
    help->show();
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

    // send resizeEvent only on show.
    CGAL::Random_points_in_square_2<Point> g(0.5);
    for(int count=0; count<25; count++) 
    {
      Point p1(*g++), p2(*g++);
      NT scale(2);
      Segment s( Point(p1.x()*scale,p1.y()*scale)  ,
                 Point(p2.x()*scale,p2.y()*scale) );
      list_of_segments.push_back(s);
      arr.insert(s);
    }

    pl_valid = false;

    something_changed();
  }
	

private:
  CGAL::Qt_widget       *widget;
  CGAL::Qt_widget_standard_toolbar
                        *stoolbar;
  Tools_toolbar         *newtoolbar;
  int                   old_state;
  Qt_layer_show_ch      testlayer;
};

#include "demo1.moc"


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
