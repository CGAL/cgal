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

#include "cgal_types2.h"

#include <fstream>
#include <stack>
#include <set>
#include <string>
#include <list>


#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_Polygon_2.h>
#include "Qt_widget_toolbar2.h"
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

typedef std::list<Cgal_Polygon> Polygon_list;
typedef std::list<Segment> Segment_list;

//global flags and variables
int                 current_state;
Polygon_list        list_of_polygons;
Arr                 arr;

Polygon_list        list_of_covering;
Segment_list        list_of_segments;

// function declaration
bool polygons_covering(Polygon_list &in_poly_list,
  Polygon_list &out_poly_list, bool intersection);


class Qt_layer_show_ch : public CGAL::Qt_widget_layer
{
public:
	
  Qt_layer_show_ch(){};


  void draw()
  {
    widget->lock();

    *widget << CGAL::GREEN;
    *widget << CGAL::LineWidth(1);
    std::list<Cgal_Polygon>::iterator itp = list_of_polygons.begin();
    while(itp!=list_of_polygons.end()){
      Cgal_Polygon::Edge_const_iterator eci = (*itp).edges_begin();
      while(eci != (*itp).edges_end() )
      {
	*widget << (*eci++);
      }
      itp++;
    }


    if( !list_of_covering.empty() )
    {
      *widget << CGAL::FillColor(CGAL::PURPLE);
      std::list<Cgal_Polygon>::iterator itp = list_of_covering.begin();
      while(itp!=list_of_covering.end())
      {
	*widget << (*itp++);
      }
    }

    *widget << CGAL::YELLOW;
    *widget << CGAL::LineWidth(3);
    std::list<Segment>::iterator its = list_of_segments.begin();
    while(its != list_of_segments.end())
    {
      *widget << (*its++);
    }

    *widget << CGAL::LineWidth(1);

    widget->unlock();
  };	


  void mousePressEvent(QMouseEvent *e)
  {
    if(e->button() == Qt::MidButton)
    {
      if( list_of_polygons.empty() )
	return;

      NT x=static_cast<NT>(widget->x_real(e->x()));
      NT y=static_cast<NT>(widget->y_real(e->y()));

      Point p(x,y);
      NT min_dist=100000000;
      Polygon_list::iterator it_closest=NULL;
      Polygon_list::iterator pit = list_of_polygons.begin();

      while(pit!=list_of_polygons.end())
      {
	Cgal_Polygon::Edge_const_iterator eit = (*pit).edges_begin();
	while(eit != (*pit).edges_end())
	{
	  NT dist = CGAL::squared_distance( p, (*eit));
	  if( dist < min_dist)
	  {
	    min_dist = dist;
	    it_closest = pit;
	  }
	  eit++;
	}
	pit++;
      }
      
      
      Arr::Curve_iterator ci = arr.curve_node_begin();
      while(ci != arr.curve_node_end() )
      {
	Cgal_Polygon::Edge_const_iterator eit=(*it_closest).edges_begin();
	while( eit!= (*it_closest).edges_end() )
	{
	  if( (*ci).curve() == (*eit) )
	  {
	    arr.remove_curve( ci );
	    break;
	  }
	  eit++;
	}
	ci++;
      }

      list_of_polygons.erase( it_closest );
      list_of_covering.clear();
      list_of_segments.clear();
      
      (*widget).redraw();
    }
  }
  
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
    QPopupMenu * algo = new QPopupMenu( this );
    menuBar()->insertItem( "&Algorithms", algo );
    algo->insertItem("&Intersection", this,
				SLOT(find_intersection()), CTRL+Key_I );
    algo->insertItem("&Union", this,
				SLOT(find_union()), CTRL+Key_U );

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
    newtoolbar = new Tools_toolbar(widget, this, &list_of_polygons);	
  
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
    list_of_polygons.clear();
    stoolbar->clear_history();
    widget->set_window(-1.1, 1.1, -1.1, 1.1); 
			// set the Visible Area to the Interval
    widget->unlock();
    something_changed();
  }

private slots:
  void get_new_object(CGAL::Object obj)
  {
    list_of_covering.clear();
    list_of_segments.clear();

    Cgal_Polygon pol;
    if(CGAL::assign(pol,obj))
    {
      if(!pol.is_simple())
      {
	QMessageBox::about( this, my_title_string,
	  "Only simple polygons are allowed.");
	return;
      }

      if(pol.orientation() == CGAL::CLOCKWISE)
	pol.reverse_orientation();

      list_of_polygons.push_back(pol);

      Cgal_Polygon::Edge_const_iterator eci = pol.edges_begin();
      while(eci != pol.edges_end() )
      {
	arr.insert( *eci++ );
      }
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

  void find_intersection()
  {
    if( list_of_polygons.empty() )
    {
      QMessageBox::about( this, my_title_string,
	"Enter some polygons.");
      return;
    }

    polygons_covering(list_of_polygons, list_of_covering, true);

    something_changed();
  }
	
  void find_union()
  {
    if( list_of_polygons.empty() )
    {
      QMessageBox::about( this, my_title_string,
	"Enter some polygons.");
      return;
    }

    polygons_covering(list_of_polygons, list_of_covering, false);

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

#include "demo2.moc"


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


//generalized face_diff function, to acount for overlaps.
int face_diff (Arr::Ccb_halfedge_circulator circ) {
  Traits t;
  int diff = 0;
  Arr::Overlap_circulator oc = circ->overlap_edges();
  do {
    if (circ->source()->point() == t.curve_source(oc->x_curve()) ) 
    diff--;     //we're inside, going outside
  else
    diff++;
  } while (++oc != circ->overlap_edges());

  return diff;
} 

// covering_DFS will compute for each face in how many polygons it is.
// It is a recursive DFS function and will be called with the unbounded 
// face after its counter has been initialized to 0.
void covering_DFS(Arr::Face_handle f) {
  Arr::Ccb_halfedge_circulator start,circ;

  // Do a recursive step for all neighbours, if any exists.
  if (f->does_outer_ccb_exist()) {
    start = circ = f->outer_ccb();
    do {
      if (circ->twin()->face()->counter < 0) {
        int diff = face_diff(circ);
        circ->twin()->face()->counter = (f->counter + diff); 
        covering_DFS(circ->twin()->face());
      }
    } while (++circ != start);
  }

  // Do a recursive step for all holes, if any exists.
  Arr::Holes_iterator hit = f->holes_begin();
  for (; Arr::Holes_iterator(hit)!=Arr::Holes_iterator(f->holes_end());
       ++hit)
  {
    start = circ = (*hit);
      do {
        if (circ->twin()->face()->counter < 0) {
          int diff = face_diff(circ);
          circ->twin()->face()->counter = (f->counter + diff); 
          covering_DFS(circ->twin()->face());        
        }
      } while (++circ != start);
  }
} 

// Convert faces of the arrangement that are in the intersection
// to polygons.
void polygons_from_faces(Arr& arr,
			 std::list<Arr::Face_iterator>& face_it_list,
			 std::list<Cgal_Polygon>& poly_list)
{
  std::list<Arr::Face_iterator>::iterator  lit;
  Cgal_Polygon                        poly;
  
  for (lit = face_it_list.begin(); lit != face_it_list.end(); lit++) {

    poly.erase(poly.vertices_begin(), poly.vertices_end());
    Arr::Ccb_halfedge_circulator cc=(*lit)->outer_ccb();
    do {
      poly.push_back(cc->curve().source());
      cc++;
    } while (cc != (*lit)->outer_ccb());
    poly_list.push_back(poly);
  }
}

// performs the extraction of data out of the processed arrangement
// if covering = 0, will perform union
// otherwise, if there are n polygons in the arrangement and covering == n
// then will perform intersection
void get_faces_with_covering(std::list<Arr::Face_iterator>& unions, 
			     int covering)
{
  Arr::Face_handle uf = arr.unbounded_face();
  uf->counter = 0;
  covering_DFS(uf);
  
  //"collecting" the union boundary faces. 
  for(Arr::Face_iterator fit = arr.faces_begin(); fit!=arr.faces_end(); ++fit) 
  {
    if (fit->counter == covering) 
    {
      unions.push_back(fit);
    }
  }
}                                                                          

void get_union()
{
  Arr::Face_handle uf = arr.unbounded_face();
  uf->counter = 0;
  covering_DFS(uf);

  Arr::Halfedge_iterator hi = arr.halfedges_begin();
  while(hi != arr.halfedges_end() )
  {
    Arr::Face_iterator fit = hi->face();
    if( fit->counter==0 )
    {
      list_of_segments.push_back( (*hi).curve() );
    }
    hi++;
  }
  
}

void clean_count()
{
  for(Arr::Face_iterator fit = arr.faces_begin(); fit!=arr.faces_end(); ++fit) 
  {
    fit->counter = -1;
  }
}

bool polygons_covering(Polygon_list &in_poly_list,
                        Polygon_list &out_poly_list, bool intersection)
{
  std::list<Arr::Face_iterator> face_it_list;

  clean_count();

  // faces with a covering two are faces that are in the intersection
  // of the two polygons.
  if(intersection)
    get_faces_with_covering(face_it_list, in_poly_list.size());
  else
  {
    get_union();
  }

  polygons_from_faces(arr, face_it_list, out_poly_list);

  if (out_poly_list.empty()) return 0; else return 1;
}

#endif // CGAL_USE_QT
