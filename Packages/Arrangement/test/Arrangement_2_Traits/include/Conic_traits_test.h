#ifndef CGAL_CONIC_TRAITS_TEST_H
#define CGAL_CONIC_TRAITS_TEST_H

#include "Base_traits_test.h"
#include <list>

static const double _error_eps = 0.001;

template <class Traits_class, class Number_type>
class Conic_traits_test : 
public Base_traits_test<Traits_class, Number_type>
{
 public:
  typedef Number_type  NT;
  typedef typename Traits_class::Point     Point;
  typedef typename Traits_class::X_curve   X_curve;
  typedef typename Traits_class::Curve     Curve;
  typedef typename Traits_class::Conic     Conic;

 public:
  
  // Constructor.
  Conic_traits_test (int argc, char** argv) :
    Base_traits_test<Traits_class, Number_type>(argc, argv) 
  {}

  // Destructor.
  virtual ~Conic_traits_test()
  {}

  // Read a curve.
  virtual void read_curve (ifstream& is, Curve& cv);
  
  // Test the make_x_monotone function.
  virtual bool make_x_monotone_wrapper (istrstream& str_line);
  
  // Test the curve_split function.
  virtual bool curve_split_wrapper (istrstream& str_line);
};

//---------------------------------------------------------------------
// Read a curve. This method is called by collect_data.
//
template <class Traits_class, class Number_type>
void Conic_traits_test<Traits_class, Number_type>::
read_curve (ifstream& is, Curve& cv)
{
  // Read a line from the input file.
  char one_line[128];

  skip_comments (is, one_line);
  istrstream str_line( one_line, 128 );

  // Get the arc type.
  char  type;
  Conic conic;

  str_line >> type;

  // An ellipse (full ellipse or a partial ellipse):
  if (type == 'f' || type == 'F' || type == 'e' || type == 'E')
  {  
    // Read the ellipse (using the format "a b x0 y0"):
    //
    //     x - x0   2      y - y0   2
    //  ( -------- )  + ( -------- )  = 1
    //       a               b
    //
    NT     a, b, x0, y0;
    
    str_line >> a >> b >> x0 >> y0;
    
    NT     a_sq = a*a;
    NT     b_sq = b*b;

    conic = Conic (b_sq, a_sq, 0,
		   -2*x0*b_sq, -2*y0*a_sq,
		   x0*x0*b_sq + y0*y0*a_sq - a_sq*b_sq);

    if (type == 'f' || type == 'F')
    {
      // Create a full ellipse.
      cv = Curve(conic);
      return;
    }
  }
  else if (type == 'h' || type == 'H')
  {
    // Read the hyperbola (using the format "a b x0 y0"):
    //
    //     x - x0   2      y - y0   2
    //  ( -------- )  - ( -------- )  = 1
    //       a               b
    //
    NT     a, b, x0, y0;
    
    str_line >> a >> b >> x0 >> y0;
    
    NT     a_sq = a*a;
    NT     b_sq = b*b;

    conic = Conic (b_sq, -a_sq, 0,
		   -2*x0*b_sq, 2*y0*a_sq,
		   x0*x0*b_sq - y0*y0*a_sq - a_sq*b_sq);  
  }
  else if (type == 'p' || type == 'P')
  {
    // Read the parabola (using the format "c x0 y0"):
    //
    //                        2
    //  4c*(y - y0) = (x - x0)
    //
    NT     c, x0, y0;
    
    str_line >> c >> x0 >> y0;
    
    conic = Conic (1, 0, 0,
		   -2*x0, -4*c,
		   x0*x0 + 4*c*y0);
  }
  else if (type == 'c' || type == 'C' || type == 'a' || type == 'A')
  {
    // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
    NT      r, s, t, u, v, w;

    str_line >> r >> s >> t >> u >> v >> w;

    conic = Conic (r, s, t, u, v, w);

    if (type == 'c' || type == 'C')
    {
      // Create a full conic (should work only for ellipses).
      cv = Curve(conic);
      return;
    }
  }
  else if (type == 's' || type == 'S')
  {
    // Read a segment, given by its endpoints (x1,y1) and (x2,y2);
    NT      x1, y1, x2, y2;

    str_line >> x1 >> y1 >> x2 >> y2;

    conic = Conic (0, 0, 0,
		   y1 - y2, x2 - x1,
		   x1*y2 - x2*y1);

    // Create the segment.
    cv = Curve(conic, Point(x1,y1), Point(x2,y2));
    return;
  }

  else
  {
    cerr << "Illegal conic type specification: " << type << "." << endl;
    return;
  }

  // Read the end points of the arc and create it.
  NT    x1, y1, x2, y2;

  str_line >> x1 >> y1 >> x2 >> y2;

  Point source (x1, y1);
  Point target (x2, y2);

  // Create the circular arc.
  cv = Curve (conic, 
	      source, target,
	      _error_eps);           // It's OK to move the endpoints a bit.
 
  return;
}

//---------------------------------------------------------------------
// Test the make_x_monotone function.
// Input case: make_x_monotone n1 n2, where: 
//     n1 - curve index in all_curves_vec
//     n2 - number of expected X-monotonian subcurves
//
template<class Traits_class, class Number_type>
bool Conic_traits_test<Traits_class, Number_type>::
make_x_monotone_wrapper (istrstream& str_line)
{
  // Read the inputs.
  int        cv_index;
  int        n_exp_curves;

  str_line >> cv_index >> n_exp_curves;

  cout << "Test: make_x_monotone( Curve" << cv_index << " ) ? " << endl;

  // Make x-monotone !
  list<X_curve> x_curves;

  tr.make_x_monotone (all_curves_vec[cv_index], x_curves);

  int           n_act_curves = static_cast<int>(x_curves.size()); 
  if (n_act_curves == n_exp_curves)
  {
    cout << "Was successful" << endl;
    return (true);
  }
  else
  {
    cout << "Was NOT successful" << endl;
    cout << "Expected " << n_exp_curves << " x-monotone curves, "
	 << "recieved " << n_act_curves << endl;
    return (false);
  }
}

//---------------------------------------------------------------------
// Test the curve_split function.
// Input case: curve_split n1 n2, where: 
//     n1 - curve index in all_curves_vec
//     n2 - point index in all_points_vec
// Does NOT take any expected result.
//
template<class Traits_class, class Number_type>
bool Conic_traits_test<Traits_class, Number_type>::
curve_split_wrapper (istrstream& str_line)
{
  // Read the inputs.
  int     cv_index, pt_index;

  str_line >> cv_index >> pt_index;
  cout << "Test: curve_split( Curve" << cv_index 
       << ", " << all_points_vec[pt_index]
       << " ) ? " << endl;

  // Check some preconditions.
  if (! tr.is_x_monotone(all_curves_vec[cv_index]))
  {
    cout << "Was NOT successful" << endl;
    cout << "Precondition fault: Input curve is not x-monotone" << endl;
    return (false);
  }

  if (tr.curve_get_point_status (all_curves_vec[cv_index], 
                                 all_points_vec[pt_index]) != tr.ON_CURVE)
  {
    cout << "Was NOT successful" << endl;
    cout << "Precondition fault: Split point is not on the curve" << endl;
    return (false);
  }

  if (all_curves_vec[cv_index].source() == all_points_vec[pt_index] ||
      all_curves_vec[cv_index].target() == all_points_vec[pt_index])
  {
    cout << "Was NOT successful" << endl;
    cout << "Precondition fault: Split point is an end point of the curve" 
	 << endl;
    return (false);
  }

  // Split the curve.
  X_curve cv1, cv2;

  tr.curve_split (all_curves_vec[cv_index], 
		  cv1, cv2, 
		  all_points_vec[pt_index]);

  // Check the results.
  if (cv1.source() != all_curves_vec[cv_index].source())
  {
    cout << "Was NOT successful" << endl;
    cout << "Source points of first parts are different" << endl;
    return (false);
  }

  if (cv1.target() != all_points_vec[pt_index])
  {
    cout << "Was NOT successful" << endl;
    cout << "Target points of first parts are different" << endl;
    return (false);
  }

  if (cv2.source() != all_points_vec[pt_index])
  {
    cout << "Was NOT successful" << endl;
    cout << "Source points of second parts are different" << endl;
    return (false);
  }

  if( cv2.target() != all_curves_vec[cv_index].target())
  {
    cout << "Was NOT successful" << endl;
    cout << "Target points of second are different" << endl;
    return (false);
  }

  cout << "Was successfull" << endl;
  return (true);
}

#endif
