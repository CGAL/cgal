#include <CGAL/basic.h>

#include <stdio.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>
#include <strstream.h>
#include <assert.h>

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
#include <CGAL/Alpha_shape_short_names_2.h>

//#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Alpha_shape_euclidean_traits_2.h>

#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>

#define CGAL::ALPHA_WINDOW_STREAM 
#include <CGAL/Alpha_shape_2.h>



#include "parse.h"



//typedef leda_integer  coord_type;
typedef double coord_type;
//typedef leda_real coord_type;
//typedef CGAL::Fixed coord_type;

typedef CGAL::Cartesian<coord_type>  Rep;
//typedef CGAL::Homogeneous<coord_type>  Rep;

//typedef CGAL::Point_2<Rep>  Point;
//typedef CGAL::Segment_2<Rep>  Segment;
typedef CGAL::Ray_2<Rep>  Ray;
typedef CGAL::Line_2<Rep>  Line;
//typedef CGAL::Triangle_2<Rep>  Triangle;



typedef CGAL::Alpha_shape_euclidean_traits_2<Rep> Gt;
typedef CGAL::Alpha_shape_euclidean_traits_2<Rep>::Point Point;
typedef CGAL::Alpha_shape_euclidean_traits_2<Rep>::Segment Segment;
typedef CGAL::Alpha_shape_euclidean_traits_2<Rep>::Triangle Triangle;


typedef CGAL::Alpha_shape_vertex_base_2<Gt> Vb;
typedef CGAL::Alpha_shape_face_base_2<Gt>  Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;


typedef CGAL::Alpha_shape_2<Gt,Tds>  Alpha_shape_2;

typedef typename Alpha_shape_2::Face  Face;
typedef typename Alpha_shape_2::Vertex Vertex;
typedef typename Alpha_shape_2::Edge Edge;
typedef typename Alpha_shape_2::Face_handle  Face_handle;
typedef typename Alpha_shape_2::Vertex_handle Vertex_handle;

typedef typename Alpha_shape_2::Face_circulator  Face_circulator;
typedef typename Alpha_shape_2::Vertex_circulator  Vertex_circulator;

typedef typename Alpha_shape_2::Locate_type Locate_type;

typedef typename Alpha_shape_2::Face_iterator  Face_iterator;
typedef typename Alpha_shape_2::Vertex_iterator  Vertex_iterator;
typedef typename Alpha_shape_2::Edge_iterator  Edge_iterator;
typedef typename Alpha_shape_2::Edge_circulator  Edge_circulator;
typedef typename Alpha_shape_2::Coord_type Coord_type;
typedef typename Alpha_shape_2::Alpha_iterator Alpha_iterator;

#define CGAL::ALPHA_WINDOW_STREAM
typedef CGAL::Window_stream  Window_stream;

void
random_input(Alpha_shape_2 &A,
	     CGAL::Window_stream &W)
  // Generate random points inside the window 
  // insert them into the alpha shape
{  
 int n = 5;
 vector<Point> V;
  V.reserve(n);

  

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

    
  for(int i = 0; i<n; i++)
    { 
      int x = rand_int(xmin,xmax);
      int y = rand_int(ymin,ymax);
      
      cerr << "Point " << x << "," << y << endl;
      Point p((double)x,(double)y);
      V.push_back(p);
     }
    cerr << "Generating " << n << " random points" << endl; 
//   vector<Point> V;
//   V.reserve(3);
//   Point p((double) 2,(double) 1);
//   V.push_back(p);
//   Point p1((double) 1,(double) 1);
//   V.push_back(p1);
//   Point p2((double) 5,(double) 2);
//   V.push_back(p2);

   n = A.make_Alpha_shape(V.begin(),V.end());

  cerr << "Inserted " << n  << " points" << endl;
}


void
container_input(Alpha_shape_2 &T,
                Window_stream &W)
{
    vector<Point> V(3);
    V[0] = Point(0, 0);
    V[1] = Point(1, 4);
    V[2] = Point(3, 4);

    //list<Vertex*> Vertices;
    int n = T.make_Alpha_shape(V.begin(), V.end());
    cerr << n << " points inserted from a vector." << endl;
    
}

void
window_input(Alpha_shape_2 &T,
             Window_stream &W)
{
    cerr << "Enter points with the left button" << endl;
    cerr << "Remove points with the middle button" << endl;
    cerr << "Right button terminates input of points" << endl;

    Point p;
    Point q(coord_type(W.xmin()-1),
            coord_type(W.ymin()-1));

    Face_handle highlight;
    Vertex_handle hv;
    vector<Point> V;
    double x, y;
    int b,n;
    while(1) {
        b = W.get_mouse(x,y);
        bool button_pressed = (b == MOUSE_BUTTON(1)) ||
                              (b == MOUSE_BUTTON(2)) ||
                              (b == MOUSE_BUTTON(3));
        p = Point(coord_type(x),
                  coord_type(y));
        bool mouse_moved = p != q;
        bool face_change = true,
             vertex_change = false;
       
        if(b == MOUSE_BUTTON(1)) {
	  W.draw_node(x,y);
           V.push_back(p);
	 
	  
            } else if(b == MOUSE_BUTTON(3)){
	  n=T.make_Alpha_shape(V.begin(),V.end());
	  cerr << n << " points inserted from a vector." << endl;
	  W.clear();
	  cerr << "Alpha optimal pour 1 cc :" << *(T.find_optimal_alpha(1)) << endl;
	  T.set_alpha(*(T.find_optimal_alpha(1)));
	  W << T;
	  cerr << "Affiche Alpha shape" << endl;
	  V.erase(V.begin(),V.end());
	  return;
        } else if(b == MOUSE_BUTTON(2))
	  {W.clear();}

        
        q = p;
    }
}


void main()
{
  Alpha_shape_2 T;
  Window_stream W(800, 800); // physical window size
// W_global = &W;

    //W.init(opt.min, opt.max, opt.min);   // logical window size
    W.init(-7,7,-7);
    W.set_show_coordinates(true);
    W << CGAL::GREEN;
    W.display();
    // random_input(T, W);
    // container_input(T,W);	
    while (1){
    window_input(T,W);}
      T.set_alpha(4);

    W.clear();
     W << T;
   cerr << "Affiche Alpha shape" << endl;
    while (1) {};

}
   
