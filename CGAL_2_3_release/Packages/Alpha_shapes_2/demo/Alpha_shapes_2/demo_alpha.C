// ======================================================================
//
// Copyright (c) 1999 The GALIA Consortium
//
// This software and related documentation is part of the
// Computational Geometry Algorithms Library (CGAL).
//
// 
// This software and documentation is provided "as-is" and without
// warranty of any kind. In no event shall the CGAL Consortium be
// liable for any damage of any kind.
// ----------------------------------------------------------------------
//
// file          : demo/Alpha_shapes_2/demo_alpha.C
// package       : Alpha_shapes_2(1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================


#include <CGAL/Cartesian.h>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <strstream>
#include <vector>
#include <list>

#include <CGAL/triangulation_assertions.h>

#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Random.h>
#define CGAL_ALPHA_WINDOW_STREAM

#include <CGAL/Alpha_shape_euclidean_traits_2.h>

#include <CGAL/Alpha_shape_vertex_base_2.h>

#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>

#include <CGAL/Real_timer.h>
#include "Parse.C"



typedef double coord_type;

// typedef leda_integer  coord_type;
// typedef leda_real coord_type;
// typedef CGAL::Fixed coord_type;

typedef CGAL::Cartesian<coord_type>  K;
// typedef CGAL::Homogeneous<coord_type>  K;

typedef K::Point_2    Point;
typedef K::Segment_2  Segment;
typedef K::Ray_2      Ray;
typedef K::Line_2     Line;
typedef K::Triangle_2 Triangle;

typedef CGAL::Alpha_shape_euclidean_traits_2<K> Gt;

typedef CGAL::Alpha_shape_vertex_base_2<Gt> Vb;

typedef CGAL::Triangulation_face_base_2<Gt> Df;
typedef CGAL::Alpha_shape_face_base_2<Gt, Df>  Fb;

typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<Gt,Tds> Triangulation_2;

typedef CGAL::Alpha_shape_2<Triangulation_2>  Alpha_shape;

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

typedef CGAL::Window_stream  Window_stream;


#ifdef CGAL_USE_CGAL_WINDOW
typedef CGAL::panel Panel;
typedef std::string String;
#else
typedef leda_panel Panel;
typedef leda_string String;
#endif



//---------------- global variables -----------------------------------

Alpha_shape* pA;

// SUN-PRO compiler problem if the object itself therefore pointer

const CGAL::Color CONTOUR_COLOR = CGAL::BLACK;
const CGAL::Color DELAUNAY_COLOR = CGAL::GREEN;
const CGAL::Color REGULARIZED_COLOR = CGAL::RED;
const CGAL::Color GENERAL_COLOR = CGAL::VIOLET;
const CGAL::Color VERTEX_COLOR = CGAL::BLUE;

const int ALPHA_MAX = 100;
const int ALPHA_MIN = 0;

CGAL::Real_timer t1;

//------------------ visualization -------------------------------------

void
any_button(Window_stream &W)
{
    double x, y;
    std::cout << "Press any button to continue" << std::endl;
    W.read_mouse(x,y);
}

//-------------------------------------------------------------------

void 
draw_vertices(const Alpha_shape& A,
		   Window_stream &W)
{
  // visualize the finite vertices
  Vertex_iterator vertex_it;

  for( vertex_it = A.vertices_begin(); 
       vertex_it != A.vertices_end();
       ++vertex_it)
    {
      Point p = vertex_it->point();
      W.draw_filled_node(p.x(), p.y());
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
      xmin = std::min( xmin, (*it).x());
      xmax = std::max( xmax, (*it).x());
      xmax = std::max( xmax,  (*it).y());
      ymin = std::min( ymin, (*it).y());
    }
 
  xmin -= 0.05*(xmax-xmin);
  ymin -= 0.05*(xmax-ymin);
  xmax += 0.25*(xmax-xmin);  // because of the panel
}

//---------------------------------------------------------------------

void
random_input(Alpha_shape &A,
	     std::vector<Point> &V, 
	     Window_stream &W,
             const Options& opt)
  // Generate random points inside the window 
  // insert them into the alpha shape
{
  int n = opt.number_of_points;
  V.reserve(n);

  std::cout << "Generating " << n << " random points" << std::endl;

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

  CGAL::Random rand;
  for(int i = 0; i<n; i++)
    { 
      int x = rand.get_int(xmin,xmax);
      int y = rand.get_int(ymin,ymax);
      Point p((double)x,(double)y);
      V.push_back(p);
     }
  t1.start();
  n = A.make_alpha_shape(V.begin(), V.end());
  t1.stop();
  std::cout << "Inserted " << n << " points in " 
	    << t1.time() << " secondes." << std::endl;
  t1.reset();
}

//---------------------------------------------------------------------

void
window_input(Alpha_shape &A,
	     std::vector<Point> &V, 
             Window_stream &W)
{
  std::cout << "Enter points with the left button" << std::endl;
  std::cout << "Right button terminates input of points" << std::endl;

  Point p;
  int n =0;

  W << VERTEX_COLOR;
  while(1) 
    {
      double x, y;
      int b = W.get_mouse(x,y);
      p = Point(Coord_type(x),
		Coord_type(y));
      
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
  std::cout << "You have entered " << V.size() << " points." << std::endl;
  t1.start();
  n = A.make_alpha_shape(V.begin(), V.end());
  t1.stop();
  std::cout << "Inserted " << n << " points in " 
	    << t1.time() << " secondes." << std::endl;
  t1.reset();

}

//--------------------------------------------------------------------

bool
file_input(std::vector<Point>& V,
	   const Options& opt)
{

  std::ifstream is(opt.finname, std::ios::in);

  if(is.fail())
    {
      std::cerr << "unable to open " << opt.finname << " for input" << std::endl;
      return false;
    }

  CGAL::set_ascii_mode(is);

  int n;
  is >> n;
  std::cout << "Reading " << n << " points" << std::endl;
  V.reserve(n);
  Point p;
  for( ; n>0 ; n--)
    {
      is >> p;
      V.push_back(p);
    }
  std::cout << "Points read" << std::endl;
  return true;
}
    

//----------------------------------------------------------------

void
file_output(const std::vector<Point>& V,
           const Options& opt)
  // the points are written in the same order as they where obtained;
  // if we woild use a vertex_iterator, this would not be the case.
{
    if(! opt.file_output)
      {
        return;
      }

    std::ofstream os(opt.foutname);
    CGAL::set_ascii_mode(os);
    
    int n = V.size();
    os << n << std::endl;

    std::vector<Point>::const_iterator it;
    for (it = V.begin(); it != V.end(); ++it)
      os << *it << std::endl;

    std::cout << n << " points written" << std::endl;
}

//-------------------------------------------------------------------

void clear_all(Alpha_shape& A,
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

void redraw(Window_stream* wp, double x0, double y0, double x1, double y1) 
{ 
  wp->flush_buffer(x0,y0,x1,y1);
}


//------------------ main -------------------------------------------

int main(int argc,  char* argv[])
{
  std::vector<Point> V;
  Alpha_shape A;
  pA = &A;
  int nn;

  //---------------- options ---------------------

  Options opt;
  parse(argc, argv, opt);
  double xmin = opt.min,  // logical window size
    xmax = opt.min,
    ymin = opt.min;

  //---------------- initialize the visualization ---------------------

  Window_stream W(opt.winx, opt.winy, "Alpha Shapes"); 
  // physical window size
 
  // buttons
  W.button("Open", 1, "Opens an input panel.");
  W.button("Save", 2, "Opens an output panel.");
  W.button("Action", 5, "Compute optimal approximation.");
  W.button("Clear", 3, "Clears point set and window.");
  W.button("Exit", 0, "Exits the program.");
  W.button("Help", 4, "Help");

  W.bool_item("Delaunay", opt.Delaunay);
  W.bool_item("Contour", opt.contour);
  W.bool_item("Regularized", opt.regularized);

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
	    Panel Pin;
	    String finname(opt.finname);


	    Pin.text_item("\\bf Get input points");
	    Pin.text_item("");
	    Pin.choice_item("", input_choice ,"File","Mouse","Random");
	    Pin.button("OK", 0);
	    Pin.button("Cancel", 1);

	    if (Pin.open(W) == 0)
	      { 
		Panel Pfin;
		
		opt.file_input = (input_choice == 0);
		switch (input_choice) 
		  {

		  case 0: 	 
		          // Get the file name
		           Pfin.string_item("Find file :", finname); //.c_str()
			  Pfin.button("OK", 0);
			  Pfin.button("Cancel", 1);
			  if (Pfin.open(W) == 1)
			    break;
#if defined(CGAL_USE_CGAL_WINDOW)			    
			  CGAL_CLIB_STD::strcpy(opt.finname, finname.c_str());  
#else			    
			  CGAL_CLIB_STD::strcpy(opt.finname, finname);
#endif

			  std::cout << opt.finname << std::endl;
	
			  clear_all(A, V, W);
			  if (!file_input(V, opt))
			    break;
				
			  xmin = xmax = ymin = opt.min;
			  get_logical_size(V.begin(), V.end(), xmin, xmax, ymin);
		
			  W.init(xmin, xmax, ymin);
			  W << VERTEX_COLOR; 
			  {
			    std::vector<Point>::iterator it;
			    for (it = V.begin(); it != V.end(); ++it)
			      W << *it;
			  }
			  t1.start();
			  nn = A.make_alpha_shape(V.begin(), V.end());
			  t1.stop();
			  std::cout << "Inserted " << nn << " points in " 
				    << t1.time() << " secondes." << std::endl;
			  t1.reset();
			  set_alpha(alpha_index);
			  W.clear();
			  W.init(xmin, xmax, ymin);
		    break;
		  case 1: clear_all(A, V, W);
		          W.init(opt.min, opt.max, opt.min);
		          window_input(A, V, W);
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
	    Panel Pout;
	    // panel P;
	    Pout.text_item("\\bf Save points");
	    Pout.text_item("");
	    String foutname(opt.foutname);
	    // Get the file name
	    Pout.string_item("Find file :", foutname); // .c_str()
	    Pout.button("OK",0);
	    Pout.button("Cancel",1);

	    if (Pout.open(W) == 0)
	      { 
		opt.file_output = true;
#if defined(CGAL_USE_CGAL_WINDOW)
		CGAL_CLIB_STD::strcpy(opt.foutname, foutname.c_str());
#else
		CGAL_CLIB_STD::strcpy(opt.foutname, foutname);
#endif
		
		file_output(V, opt);
	      }
	    break;
	  }

	case 4:
	  {
	    // help infos
	    Panel Pout;
	    Pout.text_item("Open : points input from file or mouse.");
	    Pout.text_item("");
	    Pout.text_item("Save : save points.");
	    Pout.text_item("");
	    Pout.text_item("Action : compute an optimal alpha under conditions.");
	    Pout.text_item("");
	    Pout.text_item("Clear : re-init the demo.");
	    Pout.text_item("");
	    Pout.text_item("Exit : quit the demo.");
	    Pout.text_item("");
	    Pout.button("OK",0);
	    Pout.open(W);
	    break;
	  }

	case 5:
	  {
	    // compute an optimal approximation
	    int nb_comp =1, opt_alpha_index;
	    Panel Popt;
	    // panel P;
	    Popt.text_item("\\bf Compute optimal approximation");
	    Popt.text_item("");
	    // Get number of solid components
	    Popt.int_item("Nb solid components :", nb_comp);
	    Popt.button("OK",0);
	    Popt.button("Cancel",1);

	    if ((Popt.open(W) == 0)&&(A.number_of_vertices() > 0))
	      { 
		Alpha_iterator opt_alpha_it = 
		  A.find_optimal_alpha(nb_comp);
		opt_alpha_index = int(opt_alpha_it - A.alpha_begin());
		alpha_index = (opt_alpha_index * ALPHA_MAX)/ A.number_of_alphas();
#ifdef DEBUG
		std::cout << "alpha_index :" << alpha_index << std::endl;
#endif // DEBUG
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
#ifdef DEBUG
      // write in Ascii file to test
  
      ofstream os(opt.foutname);
      CGAL::set_ascii_mode(os);
      os << A;
      std::cout << "file written" << std::endl;

      list< Vertex* > L_V;
      A.get_Alpha_shape_vertices(back_inserter(L_V));
      list<Vertex*>::iterator vertex_it = L_V.begin();
      for (; vertex_it != L_V.end(); ++vertex_it)
	std::cout << (*vertex_it)->point() << std::endl;
 
      list< pair<Face*, int> > L_E;
      A.get_Alpha_shape_edges(back_inserter(L_E));
   
#endif // DEBUG   
     */
      W.start_buffering();
      W.clear();

      if(opt.Delaunay)
	{
	  W << DELAUNAY_COLOR;
	  W << ((const Triangulation_2&) A);
	}
       
      if(opt.contour && !V.empty())
	{
	  W << CONTOUR_COLOR;
	  draw_input_contour(V.begin(), V.end() , W);
	}
      
      W << VERTEX_COLOR;
      draw_vertices(A, W);
       
      if (opt.regularized)
	{
	  A.set_mode(Alpha_shape::REGULARIZED);
	  W << REGULARIZED_COLOR;
	}
      else
	{
	  A.set_mode(Alpha_shape::GENERAL);
	  W << GENERAL_COLOR;
	}

      std::cout << "alpha_value :" << A.get_alpha() << std::endl;
      std::cout << "display Alpha Shape for alpha_index = " << alpha_index << std::endl;
      W << A;
     
      W.flush_buffer();
      W.stop_buffering();

#ifdef DEBUG
      std::cout << "nb components :" << A.number_solid_components() << std::endl;
#endif // DEBUG
    }
  return 0;
}
