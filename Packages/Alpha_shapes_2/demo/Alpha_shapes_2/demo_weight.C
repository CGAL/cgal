// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.0-I-20 $
// release_date  : $CGAL_Date: 1999/06/02 $
//
// file          : demo/Alpha_shapes_2/demo_weight.C
// package       : Alpha_shapes_2(1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

//#define CGAL_MYTRAITS

#include <CGAL/basic.h>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <strstream>

// #include <cassert>
#include <CGAL/triangulation_assertions.h>

// #include <CGAL/Cartesian.h>
// //#include <CGAL/Homogeneous.h>
// //#include <CGAL/Integer.h>
// //#include <CGAL/Rational.h>
// //#include <CGAL/Fixed.h>
// //#include <CGAL/Real.h>
// #include <CGAL/squared_distance_2.h>   // to avoid a g++ problem
// #include <CGAL/Point_2.h>
// #include <CGAL/predicates_on_points_2.h>
// #include <CGAL/Triangle_2.h>
// #include <CGAL/Segment_2.h>

//#include <CGAL/Weighted_alpha_shape_short_names_2.h>
//#include <CGAL/Triangulation_short_names_2.h>

#include <CGAL/Alpha_shape_vertex_base_2.h>

#include <list>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Weighted_alpha_shape_face_base_2.h>

//#include <CGAL/IO/Window_stream.h>
#define CGAL_WEIGHTED_ALPHA_WINDOW_STREAM

#ifndef CGAL_MYTRAITS
#include <CGAL/Weighted_alpha_shape_euclidean_traits_2.h>
#else
#include <CGAL/Alpha_shape_euclidean_mytraits_2.h>
#endif

#include <CGAL/Circle_2.h>
#include <CGAL/Weighted_point.h>
#include <CGAL/Weighted_alpha_shape_2.h>

#include "Parse_weight.C"
#include "Timing.C"


//typedef leda_integer  coord_type;
//typedef double coord_type;
//typedef leda_real coord_type;
//typedef CGAL::Fixed coord_type;

typedef CGAL::Cartesian<double>  Rep;
//typedef CGAL::Homogeneous<coord_type>  Rep;

typedef CGAL::Point_2<Rep> Point_base;
typedef CGAL::Weighted_point<Point_base,double>  Point;
typedef CGAL::Segment_2<Rep>  Segment;
typedef CGAL::Ray_2<Rep>  Ray;
typedef CGAL::Line_2<Rep>  Line;
typedef CGAL::Triangle_2<Rep>  Triangle;

typedef CGAL::Weighted_alpha_shape_euclidean_traits_2<Rep> Gt;



typedef CGAL::Alpha_shape_vertex_base_2<Gt> Vb;
typedef CGAL::Weighted_alpha_shape_face_base_2<Gt>  Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
typedef CGAL::Regular_triangulation_2<Gt,Tds> Dtriangulation_2;

typedef CGAL::Weighted_alpha_shape_2<Gt,Tds>  Alpha_shape_2;

typedef Alpha_shape_2::Face  Face;
typedef Alpha_shape_2::Vertex Vertex;
typedef Alpha_shape_2::Edge Edge;
typedef Alpha_shape_2::Face_handle  Face_handle;
typedef Alpha_shape_2::Vertex_handle Vertex_handle;

typedef Alpha_shape_2::Face_circulator  Face_circulator;
typedef Alpha_shape_2::Vertex_circulator  Vertex_circulator;

typedef Alpha_shape_2::Locate_type Locate_type;

typedef Alpha_shape_2::Face_iterator  Face_iterator;
typedef Alpha_shape_2::Vertex_iterator  Vertex_iterator;
typedef Alpha_shape_2::Edge_iterator  Edge_iterator;
typedef Alpha_shape_2::Edge_circulator  Edge_circulator;
//typedef Alpha_shape_2::Line_face_circulator  Line_face_circulator;
//typedef Alpha_shape_2::Coord_type Coord_type;
typedef Alpha_shape_2::Alpha_iterator Alpha_iterator;

typedef CGAL::Window_stream  Window_stream;
//---------------- global variables -----------------------------------

Alpha_shape_2* pA;

// SUN-PRO compiler problem if the object itself therefore pointer

const CGAL::Color CONTOUR_COLOR = CGAL::BLACK;
const CGAL::Color DELAUNAY_COLOR = CGAL::GREEN;
const CGAL::Color REGULARIZED_COLOR = CGAL::RED;
const CGAL::Color GENERAL_COLOR = CGAL::VIOLET;
const CGAL::Color VERTEX_COLOR = CGAL::BLUE;

const int ALPHA_MAX = 100;
const int ALPHA_MIN = 0;

//------------------ visualization -------------------------------------

void
any_button(Window_stream &W)
{
    double x, y;
    std::cerr << "Press any button to continue" << std::endl;
    W.read_mouse(x,y);
}

//-------------------------------------------------------------------

void 
draw_vertices(const Alpha_shape_2& A,
		   Window_stream &W, bool option)
{
  // visualize the finite vertices
  Vertex_iterator vertex_it;

  for( vertex_it = A.vertices_begin(); 
       vertex_it != A.vertices_end();
       ++vertex_it)
    {
      Point p = vertex_it->point();
      if (option)
	{ W << CGAL::Circle_2<Rep>(p.point(),max(p.weight(),DBL_MIN)); }
      else
	{ W.draw_filled_node(p.x(), p.y()); }
    }
}


//-------------------------------------------------------------------

template<class InputIterator>
void
draw_input_contour(InputIterator first,  
		   InputIterator last,
		   Window_stream &W)
{
  // visualize the contours

  InputIterator it;
  Point pred = *(first);

  for( it = first; 
       it != last;
       ++it)
    {
      W << Segment(pred, *it);
      pred = *it;
    }
}

//-------------------------------------------------------------------

template<class InputIterator>
void
get_logical_size(InputIterator first,  
		 InputIterator last,
		 double& xmin,
		 double& xmax,
		 double& ymin)
{
  // visualize the contours

  InputIterator it;
 
  for( it = first; 
       it != last;
       ++it)
    {
      xmin = min( xmin, (*it).x());
      xmax = max( xmax, (*it).x());
      xmax = max( xmax,  (*it).y());
      ymin = min( ymin, (*it).y());
    }
 
  xmin -= 0.05*(xmax-xmin);
  ymin -= 0.05*(xmax-ymin);
  xmax += 0.35*(xmax-xmin);  // because of the panel
}

//---------------------------------------------------------------------

void
random_input(Alpha_shape_2 &A,
	     std::vector<Point> &VV, 
	     Window_stream &W,
             const Options& opt)
  // Generate random points inside the window 
  // insert them into the alpha shape
{
  std::vector<Point> V;
  int n = opt.number_of_points;
  V.reserve(n);

  std::cerr << "Generating " << n << " random points" << std::endl;

  double x0 = W.xmin();
  double y0 = W.ymin();
  double x1 = W.xmax();
  double y1 = W.ymax();

  double dx = x1 - x0;
  double dy = y1 - y0;
 
  int xmin = int(x0 + 0.1*dx);
  int xmax = int(x1 - 0.1*dx);
  int ymin = int(y0 + 0.1*dy);
  int ymax = int(y1 - 0.1*dy);

  W << VERTEX_COLOR;

  for(int i = 0; i<n; i++)
    { 
      int x = rand_int(xmin,xmax);
      int y = rand_int(ymin,ymax);
      Point p(Point::Point((double)x,(double)y));
      V.push_back(p);
     }
  start_timing();
  if (opt.init)
    { VV=A.initialize_weighted_points_to_the_nearest_voronoi_edge(V.begin(), V.end()); }
  else
    { VV=A.initialize_weighted_points_to_the_nearest_vertex(V.begin(), V.end()); }
  n = A.make_Alpha_shape(VV.begin(), VV.end());
  end_timing(1);
  std::cerr << "Inserted " << n  << " points" << std::endl;
}

//---------------------------------------------------------------------

void
window_input(Alpha_shape_2 &A,
	     std::vector<Point> &VV, 
             Window_stream &W,
             const Options& opt)
{
  std::cerr << "Enter points with the left button" << std::endl;
  std::cerr << "Right button terminates input of points" << endl;

  std::vector<Point> V;
  Point p;
  int n =0;

  W << VERTEX_COLOR;
  while(1) 
    {
      double x, y;
      int b = W.get_mouse(x,y);
      p = Point(Point::Point(double (x),
			     double (y)));
      
      if(b == MOUSE_BUTTON(1))
	{
	  V.push_back(p);
	  W.draw_filled_node(x, y); 
        }  
      else 
	if(b == MOUSE_BUTTON(3))
	  {
	  // we are done. 
	  break;
	  }
    }
  std::cerr << "You have entered " << V.size() << " points." << endl;
  start_timing();

  if (opt.init)
    { VV=A.initialize_weighted_points_to_the_nearest_voronoi_edge(V.begin(), V.end()); }
  else
    { VV=A.initialize_weighted_points_to_the_nearest_vertex(V.begin(), V.end()); }
  n = A.make_Alpha_shape(VV.begin(), VV.end());
  end_timing(1);
  std::cerr << "Inserted " << n  << " points" << endl;

}



//--------------------------------------------------------------------

bool
file_input(Alpha_shape_2& A, 
	   std::vector<Point>& VV,
	   Window_stream &W,
           const Options& opt)
{
  std::vector<Point> V;
  ifstream is(opt.finname, ios::in, filebuf::openprot);

  if(is.fail())
    {
      std::cerr << "unable to open " << opt.finname << " for input" << endl;
      return false;
    }

  CGAL::set_ascii_mode(is);

  int n;
  is >> n;
  std::cerr << "Reading " << n << " points" << endl;
  V.reserve(n);
  Point_base p;
  for( ; n>0 ; n--)
    {
      is >> p;
      V.push_back(Point(p));
    }
  if (opt.init)
    { VV=A.initialize_weighted_points_to_the_nearest_voronoi_edge(V.begin(), V.end()); }
  else
    { VV=A.initialize_weighted_points_to_the_nearest_vertex(V.begin(), V.end()); }
  
  std::cerr << "Points read" << endl;
  return true;
}
    

//----------------------------------------------------------------

void
file_output(std::vector<Point>& V,
           const Options& opt)
  // the points are written in the same order as they where obtained;
  // if we woild use a vertex_iterator, this would not be the case.
{
  
    if(! opt.file_output)
      {
        return;
      }
    
    ofstream os(opt.foutname);
    CGAL::set_ascii_mode(os);
    
    int n = V.size();
    os << n << endl;

    std::vector<Point>::const_iterator it;
    for (it = V.begin(); it != V.end(); ++it)
      os << Point_base (*it) << endl;

    std::cerr << n << " points written" << endl;
}

//-------------------------------------------------------------------

void clear_all(Alpha_shape_2& A,
	       std::vector<Point>& V, 
	       Window_stream& W)
  // clears the window and the structures
{ 
  W.start_buffering();
  W.clear();
  W.flush_buffer();
  W.stop_buffering();

  A.clear();
  V.erase(V.begin(), V.end());
}

//-------------------------------------------------------------------


void set_alpha(int alpha_index)
{
  // alpha corresponds to an index
  if (pA->number_of_alphas() > 0)
  {
    if (alpha_index < ALPHA_MAX)
      {
	int n = (alpha_index * pA->number_of_alphas())/ ALPHA_MAX;
	pA->set_alpha(pA->get_nth_alpha(n));
      }
    else
      {
	Alpha_iterator alpha_end_it = pA->alpha_end();
	pA->set_alpha((*(--alpha_end_it))+1);
      }
  }
  else
    pA->set_alpha(ALPHA_MIN);
}

//-------------------------------------------------------------------

void redraw(leda_window* wp, double x0, double y0, double x1, double y1) 
{ 
  wp->flush_buffer(x0,y0,x1,y1);
}



//------------------ main -------------------------------------------

int main(int argc,  char* argv[])
{
  std::vector<Point> V;
  Alpha_shape_2 A;
  pA = &A;
  int nn;

  //---------------- options ---------------------

  Options opt;
  parse(argc, argv, opt);
  double xmin = opt.min,  // logical window size
    xmax = opt.min,
    ymin = opt.min;

  //---------------- initialize the visualization ---------------------

  Window_stream W(opt.winx, opt.winy, " Weighted Alpha Shapes"); 
  // physical window size
 
  // buttons
  W.button("Open", 1, "Opens an input panel.");
  W.button("Save", 2, "Opens an output panel.");
  W.button("Action", 5, "Compute optimal approximation.");
  W.button("Clear", 3, "Clears point set and window.");
  W.button("Exit", 0, "Exits the program.");
  W.button("Help", 4, "Help");

  W.bool_item("Regular", opt.Delaunay);
  W.bool_item("Contour", opt.contour);
  W.bool_item("Regularized", opt.regularized);
  W.bool_item("Show Weights", opt.weight);
  W.bool_item("Init C/ppv", opt.init);

  int alpha_index=(ALPHA_MAX-ALPHA_MIN)/2;
  W.int_item("Alpha index", alpha_index, ALPHA_MIN, ALPHA_MAX, set_alpha);

  
  W.display();
  W.init(opt.min, opt.max, opt.min);  
  // logical window size
  W.set_redraw(redraw);
  W.set_node_width(2);
  W.set_show_coordinates(true);
  
  for(;;)
    {
      double x, y;
      int but = W.read_mouse(x,y);

      if (but == 0)
	exit(0);

      switch (but) 
	{
	case 1:
	  { // get input points
	    int input_choice;
	    leda_panel Pin;
	    leda_string finname(opt.finname);

	    Pin.text_item("\\bf Get input points");
	    Pin.text_item("");
	    Pin.choice_item("", input_choice ,"File","Mouse","Random");
	    Pin.button("OK",0);
	    Pin.button("Cancel",1);

	    if (Pin.open(W) == 0)
	      { 
		leda_panel Pfin;
		
		opt.file_input = (input_choice == 0);
		switch (input_choice) 
		  {
		  case 0: 	 
		          // Get the file name
		          Pfin.string_item("Find file :", finname);
			  Pfin.button("OK",0);
			  Pfin.button("Cancel",1);
			  if (Pfin.open(W) == 1)
			    break;
			  strcpy(opt.finname, finname);
			  std::cerr << opt.finname << std::endl;
	
			  clear_all(A, V, W);
			  if (!file_input(A, V, W, opt))
			    break;
				
			  xmin = xmax = ymin = opt.min;
			  get_logical_size(V.begin(), V.end(), xmin, xmax, ymin);
		
			  W.init(xmin, xmax, ymin);
			  W << VERTEX_COLOR; 
			  std::vector<Point>::const_iterator it;
			  for (it = V.begin(); it != V.end(); ++it)
			    W << *it;
			  
			  start_timing();
			  nn = A.make_Alpha_shape(V.begin(), V.end());
			  end_timing(1);
			  std::cerr << "Inserted " << nn  << " points" << std::endl;
			  set_alpha(alpha_index);
			  W.clear();
			  W.init(xmin, xmax, ymin);
		    break;
		  case 1: clear_all(A, V, W);
		          W.init(opt.min, opt.max, opt.min);
		          window_input(A, V, W, opt);
			  set_alpha(alpha_index);
		    break;
		  case 2: clear_all(A, V, W);
		          W.init(opt.min, opt.max, opt.min);
		          random_input(A, V, W, opt);
			  set_alpha(alpha_index);
                    break;
		  default: 
		    break;
		  }

		
		
	      }

	    break;
	  }
   
	case 2:
	  {
	    // write points
	    leda_panel Pout;
	    // panel P;
	    Pout.text_item("\\bf Save points");
	    Pout.text_item("");
	    leda_string foutname(opt.foutname);
	    // Get the file name
	    Pout.string_item("Find file :", foutname);
	    Pout.button("OK",0);
	    Pout.button("Cancel",1);

	    if (Pout.open(W) == 0)
	      { 
		opt.file_output = true;
		strcpy(opt.foutname, foutname);
		file_output(V, opt);
	      }
	    break;
	  }

	case 5:
	  {
	    // compute an optimal approximation
	    int nb_comp =1, opt_alpha_index;
	    leda_panel Popt;
	    // panel P;
	    Popt.text_item("\\bf Compute optimal approximation");
	    Popt.text_item("");
	    // Get number of solid components
	    Popt.int_item("Nb solid components :", nb_comp);
	    Popt.button("OK",0);
	    Popt.button("Cancel",1);

	    if (Popt.open(W) == 0)
	      { 
		Alpha_iterator opt_alpha_it = 
		  A.find_optimal_alpha(nb_comp);
		opt_alpha_index = int(opt_alpha_it - A.alpha_begin());
		alpha_index = (opt_alpha_index * ALPHA_MAX)/ A.number_of_alphas();
#ifndef NDEBUG
		std::cerr << "alpha_index :" << alpha_index << std::endl;
#endif // NDEBUG
		A.set_alpha((*opt_alpha_it + *(opt_alpha_it+1))/2);
	      }
	    break;
	  }

	case 3: 
	  { // clear
	    clear_all(A, V, W);
	    alpha_index = (ALPHA_MAX-ALPHA_MIN)/2;
	    continue;
	  } 
	default:
	  break;
	};

      /*      // draw
#ifndef NDEBUG
      // write in Ascii file to test
  
      ofstream os(opt.foutname);
      CGAL::set_ascii_mode(os);
      os << A;
      std::cerr << "file written" << std::endl;

      std::list< Vertex* > L_V;
      A.get_Alpha_shape_vertices(back_inserter(L_V));
      std::list<Vertex*>::iterator vertex_it = L_V.begin();
      for (; vertex_it != L_V.end(); ++vertex_it)
	std::cerr << (*vertex_it)->point() << std::endl;
 
      std::list< std::pair<Face*, int> > L_E;
      A.get_Alpha_shape_edges(back_inserter(L_E));
   
#endif // NDEBUG   
     */
      
      
      W.start_buffering();
      W.clear();
      
      if (opt.weight) 
	{ 
	  W << REGULARIZED_COLOR;
	  W << CGAL::Circle_2<Rep> (Point_base ((xmax-(xmax-xmin)/15), (ymin+(xmax-xmin)/15)),max(A.get_alpha(),DBL_MIN)); 
	}
      
      if(opt.Delaunay)
	{
	  W << DELAUNAY_COLOR;
	  W << ((const Dtriangulation_2&) A);
	}
       
      if(opt.contour && !V.empty())
	{
	  W << CONTOUR_COLOR;
	  draw_input_contour(V.begin(), V.end() , W);
	}
      
      W << VERTEX_COLOR;
      draw_vertices(A, W,opt.weight);
       
      if (opt.regularized)
	{
	  A.set_mode(Alpha_shape_2::REGULARIZED);
	  W << REGULARIZED_COLOR;
	}
      else
	{
	  A.set_mode(Alpha_shape_2::GENERAL);
	  W << GENERAL_COLOR;
	}
      std::cerr << "display Weighted Alpha Shape for alpha_index = " << alpha_index << std::endl;
      W << A;

      W.flush_buffer();
      W.stop_buffering();

#ifndef NDEBUG
      std::cerr << "nb components :" << A.number_solid_components() << std::endl;
#endif // NDEBUG
    }
}
