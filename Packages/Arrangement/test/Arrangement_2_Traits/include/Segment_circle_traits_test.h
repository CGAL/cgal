#ifndef CGAL_SEGMENT_CIRCLE_TRAITS_TEST_H
#define CGAL_SEGMENT_CIRCLE_TRAITS_TEST_H

#include "Base_traits_test.h"
#include <list>

template <class Traits_class, class Number_type>
class Segment_circle_traits_test : 
public Base_traits_test<Traits_class, Number_type>
{
 public:
  typedef Number_type  NT;
  typedef typename Traits_class::Point     Point;
  typedef typename Traits_class::X_curve   X_curve;
  typedef typename Traits_class::Curve     Curve;
  typedef typename Traits_class::Segment   Segment;
  typedef typename Traits_class::Circle    Circle;

 public:
  
  // Constructor.
  Segment_circle_traits_test (int argc, char** argv) :
    Base_traits_test<Traits_class, Number_type>(argc, argv) 
  {}

  // Destructor.
  virtual ~Segment_circle_traits_test()
  {}

  // Read a curve.
  virtual void read_curve (std::ifstream& is, Curve& cv);
  
  // Test the make_x_monotone function.
  virtual bool make_x_monotone_wrapper (std::istrstream& str_line);
  
  // Test the curve_split function.
  virtual bool curve_split_wrapper (std::istrstream& str_line);
};

//---------------------------------------------------------------------
// Read a curve. This method is called by collect_data.
//
template <class Traits_class, class Number_type>
void Segment_circle_traits_test<Traits_class, Number_type>::
read_curve (std::ifstream& is, Curve& cv)
{
  // Read a line from the input file.
  char one_line[128];

  skip_comments (is, one_line);
  std::istrstream str_line( one_line, 128 );

  // Get the arc type.
  char type;

  str_line >> type;
  
  // A full circle (c) or a circular arc (a):
  if (type == 'c' || type == 'C' || type == 'a' || type == 'A')
  {  
    // Read the circle, using the format "x0 y0 r^2"
    NT     x0, y0, r2;
    
    str_line >> x0 >> y0 >> r2;
    
    Circle circle (Point (x0, y0), r2);

    if (type == 'c' || type == 'C')
    {
      // Create a full circle.
      cv = Curve(circle);  
    }
    else
    {
      // Read the end points of the circular arc.
      NT    x1, y1, x2, y2;

      str_line >> x1 >> y1 >> x2 >> y2;

      if ((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0) != r2)
	y1 = CGAL::sqrt(r2 - (x1 - x0)*(x1 - x0)) + y0;

      if ((x2 - x0)*(x2 - x0) + (y2 - y0)*(y2 - y0) != r2)
	y2 = CGAL::sqrt(r2 - (x2 - x0)*(x2 - x0)) + y0;

      Point source (x1, y1);
      Point target (x2, y2);

      // Create the circular arc.
      cv = Curve (circle, source, target);
    }
  }
  else if (type == 's' || type == 'S')
  {
    // Read the end points of the segment.
    NT    x1, y1, x2, y2;

    str_line >> x1 >> y1 >> x2 >> y2;

    Point source (x1, y1);
    Point target (x2, y2);
   
    cv = Curve (Segment (source, target));
  }
  else
  {
    // Illegal type!
    std::cout << "Failed to read the curve specified by '" 
	      << one_line << "'" << std::endl;
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
bool Segment_circle_traits_test<Traits_class, Number_type>::
make_x_monotone_wrapper (std::istrstream& str_line)
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
bool Segment_circle_traits_test<Traits_class, Number_type>::
curve_split_wrapper (std::istrstream& str_line)
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
    std::cout << "Precondition fault: Input curve is not x-monotone" << std::endl;
    return (false);
  }

  if (tr.curve_get_point_status (all_curves_vec[cv_index], 
                                 all_points_vec[pt_index]) != tr.ON_CURVE)
  {
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "Precondition fault: Split point is not on the curve" << std::endl;
    return (false);
  }

  if (all_curves_vec[cv_index].source() == all_points_vec[pt_index] ||
      all_curves_vec[cv_index].target() == all_points_vec[pt_index])
  {
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "Precondition fault: Split point is an end point of the curve" 
	 << std::endl;
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
