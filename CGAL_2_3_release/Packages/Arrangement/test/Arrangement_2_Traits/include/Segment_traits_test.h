#include "Base_traits_test.h"

template< class Traits_class, class Number_type >
class Segment_traits_test : public Base_traits_test< Traits_class, Number_type >{
public:
  typedef Number_type  NT;
  typedef typename Traits_class::Point     Point;
  typedef typename Traits_class::X_curve   X_curve;
  typedef typename Traits_class::Curve     Curve;
public:
  Segment_traits_test( int argc, char** argv );
  virtual void read_curve( std::ifstream& is, Curve& cv );
  virtual bool make_x_monotone_wrapper( std::istrstream& strLine );
  virtual bool curve_split_wrapper( std::istrstream& strLine );
  ~Segment_traits_test();
};

//------------------------------------------------------------------------------
/*
  Constructor. Nothing to do. Just calls super class ctor
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
Segment_traits_test< Traits_class, Number_type >::
Segment_traits_test( int argc, char** argv ) :
 Base_traits_test< Traits_class, Number_type >(argc, argv) {}; 
//------------------------------------------------------------------------------
/*
  Destructor. Nothing to do. Implements super class virtual dtor.
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
Segment_traits_test< Traits_class, Number_type >::
~Segment_traits_test( ){}
//------------------------------------------------------------------------------
/*
  Reads one curve. This method is called by collect_data.
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
void Segment_traits_test< Traits_class, Number_type >::
read_curve( std::ifstream& is, Curve& cv ){
  char one_line[128];
  NT x1,y1, x2, y2;

  skip_comments( is, one_line );
  std::istrstream strLine( one_line, 128 );
  strLine >> x1 >> y1; 
  skip_comments( is, one_line );
  std::istrstream strLine2( one_line, 128 );
  strLine2 >> x2 >> y2;
  cv = Curve( Point( x1,y1 ), Point( x2, y2 ) );
}
//------------------------------------------------------------------------------
/*
  input case:
  make_x_monotone n1 n2, where 
  n1 - curve index in all_curves_vec
  n2 - number of expected X-monotonian subcurves
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
bool Segment_traits_test< Traits_class, Number_type >::
make_x_monotone_wrapper( std::istrstream& strLine ){
  std::cout << "Test: make_x_monotone - nothing to do in segment case" << std::endl;
  std::cout << "Was successful" << std::endl;
  return true;
}
//------------------------------------------------------------------------------
/*
  input case:
  curve_split n1 n2, where
  n1 - curve index in all_curves_vec
  n2 - point index in all_points_vec
  Does NOT take any expected result
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
bool Segment_traits_test< Traits_class, Number_type >::
curve_split_wrapper( std::istrstream& strLine ){
  int index1, index2;
  X_curve cv1, cv2;

  strLine >> index1 >> index2;
  std::cout << "Test: curve_split( Curve" << index1 << ", " << all_points_vec[index2]
       << " ) ? " << std::endl;
  if( !tr.is_x_monotone( all_curves_vec[index1] ) ){
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "Precondition fault: Input curve is not x-monotone" << std::endl;
    return false;
  }
  if( tr.curve_get_point_status( all_curves_vec[index1], 
                                 all_points_vec[index2] ) != tr.ON_CURVE ){
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "Precondition fault: Split point is not on the curve " << std::endl;
    return false;
  }   
  if( all_curves_vec[index1].source() == all_points_vec[index2] ||
      all_curves_vec[index1].target() == all_points_vec[index2]  ){
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "Precondition fault: Split point is the end point of the curve " << std::endl;
    return false;
  }
  tr.curve_split( all_curves_vec[index1], cv1, cv2, all_points_vec[index2] );
  if( cv1.source() != all_curves_vec[index1].source() )
  {
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "Source points of first parts are different" << std::endl;
    std::cout << "The original source point: " << cv1.source() << std::endl;
    std::cout << "The obtained source point: " << all_curves_vec[index1].source() << std::endl;
    return false;
  }
  if( cv1.target() != all_points_vec[index2] )
  {
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "Target points of first parts are different" << std::endl;
    std::cout << "The original target point: " << cv1.target() << std::endl;
    std::cout << "The obtained target point: " << all_curves_vec[index1].target() << std::endl;
    return false;
  }
  if( cv2.source() != all_points_vec[index2] )
  {
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "Source points of second parts are different" << std::endl;
    std::cout << "The original source point: " << cv1.source() << std::endl;
    std::cout << "The obtained source point: " << all_curves_vec[index1].source() << std::endl;
    return false;
  }
  if( cv2.target() != all_curves_vec[index1].target() )
  {
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "Target points of second are different" << std::endl;
    std::cout << "The original target point: " << cv1.target() << std::endl;
    std::cout << "The obtained target point: " << all_curves_vec[index1].target() << std::endl;
    return false;
  }

  std::cout << "Was successfull" << std::endl;
  return true;
}
//------------------------------------------------------------------------------
