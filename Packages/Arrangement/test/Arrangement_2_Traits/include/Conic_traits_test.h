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

  typedef typename Traits_class::Point_2        Point;
  typedef typename Traits_class::Segment_2      Segment;
  typedef typename Traits_class::Circle_2       Circle;
  typedef typename Traits_class::X_curve_2      X_curve;
  typedef typename Traits_class::Curve_2        Curve;

 public:
  
  // Constructor.
  Conic_traits_test (int argc, char** argv) :
    Base_traits_test<Traits_class, Number_type>(argc, argv) 
  {}

  // Destructor.
  virtual ~Conic_traits_test()
  {}

  // Read a curve.
  virtual void read_curve (std::ifstream & is, Curve & cv);
  
  // Test the make_x_monotone function.
  virtual bool make_x_monotone_wrapper (std::istringstream& str_line);
  
  // Test the curve_split function.
  virtual bool curve_split_wrapper (std::istringstream& str_line);
};

//---------------------------------------------------------------------
// Read a curve. This method is called by collect_data.
//
template <class Traits_class, class Number_type>
void Conic_traits_test<Traits_class, Number_type>::
read_curve (std::ifstream & is, Curve & cv)
{
  // Read a line from the input file.
  char one_line[128];

  skip_comments (is, one_line);
  std::string stringvalues(one_line);
  std::istringstream str_line(stringvalues, std::istringstream::in);

  // Get the arc type.
  char     type;
  bool     is_circle = false;              // Is this a circle.
  Circle   circle;
  NT       r, s, t, u, v, w;               // The conic coefficients.

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

    if (a == b)
    {
      is_circle = true;
      circle = Circle (Point (x0, y0), a*b, CGAL::CLOCKWISE);
    }
    else
    {
      r = b_sq;
      s = a_sq;
      t = 0;
      u = -2*x0*b_sq;
      v = -2*y0*a_sq;
      w = x0*x0*b_sq + y0*y0*a_sq - a_sq*b_sq;
    }

    if (type == 'f' || type == 'F')
    {
      // Create a full ellipse (or circle).
      if (is_circle)
	cv = Curve (circle);
      else
	cv = Curve (r, s, t, u, v, w);

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

    r = b_sq;
    s= -a_sq;
    t = 0;
    u = -2*x0*b_sq;
    v = 2*y0*a_sq;
    w = x0*x0*b_sq - y0*y0*a_sq - a_sq*b_sq;  
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
    
    r = 1;
    s = 0;
    t = 0;
    u = -2*x0;
    v = -4*c;
    w = x0*x0 + 4*c*y0;
  }
  else if (type == 'c' || type == 'C' || type == 'a' || type == 'A')
  {
    // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
    str_line >> r >> s >> t >> u >> v >> w;

    if (type == 'c' || type == 'C')
    {
      // Create a full conic (should work only for ellipses).
      cv = Curve (r, s, t, u, v, w);
      return;
    }
  }
  else if (type == 's' || type == 'S')
  {
    // Read a segment, given by its endpoints (x1,y1) and (x2,y2);
    NT      x1, y1, x2, y2;

    str_line >> x1 >> y1 >> x2 >> y2;

    Point   source (x1, y1);
    Point   target (x2, y2);
    Segment segment (source, target);

    // Create the segment.
    cv = Curve (segment);
    return;
  }
  else if (type == 'i' || type == 'I')
  {
    // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
    str_line >> r >> s >> t >> u >> v >> w;

    // Read the approximated source, along with a general conic 
    // <r_1,s_1,t_1,u_1,v_1,w_1> whose intersection with <r,s,t,u,v,w>
    // defines the source.
    NT     r1, s1, t1, u1, v1, w1;
    NT     x1, y1;

    str_line >> x1 >> y1;
    str_line >> r1 >> s1 >> t1 >> u1 >> v1 >> w1;

    Point   app_source (x1, y1);

    // Read the approximated target, along with a general conic 
    // <r_2,s_2,t_2,u_2,v_2,w_2> whose intersection with <r,s,t,u,v,w>
    // defines the target.
    NT     r2, s2, t2, u2, v2, w2;
    NT     x2, y2;

    str_line >> x2 >> y2;
    str_line >> r2 >> s2 >> t2 >> u2 >> v2 >> w2;

    Point   app_target (x2, y2);

    // Create the conic arc.
    cv = Curve (r, s, t, u, v ,w,
		app_source, r1, s1, t1, u1, v1, w1,
		app_target, r2, s2, t2, u2, v2, w2);
    return;
  }
  else
  {
    std::cerr << "Illegal conic type specification: " << type << "."
              << std::endl;
    return;
  }

  // Read the end points of the arc and create it.
  NT    x1, y1, x2, y2;

  str_line >> x1 >> y1 >> x2 >> y2;

  Point source (x1, y1);
  Point target (x2, y2);

  // Create the conic (or circular) arc.
  if (is_circle)
  {
    cv = Curve (circle,
		  source, target);
  }
  else
  {
    cv = Curve (r, s, t, u, v, w,
		  source, target);
  }
 
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
make_x_monotone_wrapper (std::istringstream& str_line)
{
  // Read the inputs.
  int        cv_index;
  int        n_exp_curves;

  str_line >> cv_index >> n_exp_curves;

  std::cout << "Test: make_x_monotone( Curve" << cv_index << " ) ? " 
	    << std::endl;

  // Make x-monotone !
  std::list<X_curve> x_curves;

  tr.make_x_monotone (all_curves_vec[cv_index], x_curves);

  int           n_act_curves = static_cast<int>(x_curves.size()); 
  if (n_act_curves == n_exp_curves)
  {
    std::cout << "Was successful" << std::endl;
    return (true);
  }
  else
  {
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "Expected " << n_exp_curves << " x-monotone curves, "
	 << "recieved " << n_act_curves << std::endl;
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
curve_split_wrapper (std::istringstream& str_line)
{
  // Read the inputs.
  int     cv_index, pt_index;

  str_line >> cv_index >> pt_index;
  std::cout << "Test: curve_split( Curve" << cv_index 
       << ", " << all_points_vec[pt_index]
       << " ) ? " << std::endl;

  // Check some preconditions.
  if (! tr.is_x_monotone(all_curves_vec[cv_index]))
  {
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "Precondition fault: Input curve is not x-monotone" 
	      << std::endl;
    return (false);
  }

  if (! tr.curve_is_in_x_range (all_curves_vec[cv_index], 
				all_points_vec[pt_index]) ||
      tr.curve_get_point_status (all_curves_vec[cv_index], 
                                 all_points_vec[pt_index]) != CGAL::EQUAL)
  {
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "Precondition fault: Split point is not on the curve" 
	      << std::endl;
    return (false);
  }

  if (all_curves_vec[cv_index].source() == all_points_vec[pt_index] ||
      all_curves_vec[cv_index].target() == all_points_vec[pt_index])
  {
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "Precondition fault: "
	      << "Split point is an end point of the curve" << std::endl;
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
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "Source points of first parts are different" << std::endl;
    return (false);
  }

  if (cv1.target() != all_points_vec[pt_index])
  {
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "Target points of first parts are different" << std::endl;
    return (false);
  }

  if (cv2.source() != all_points_vec[pt_index])
  {
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "Source points of second parts are different" << std::endl;
    return (false);
  }

  if( cv2.target() != all_curves_vec[cv_index].target())
  {
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "Target points of second are different" << std::endl;
    return (false);
  }

  std::cout << "Was successfull" << std::endl;
  return (true);
}

#endif
