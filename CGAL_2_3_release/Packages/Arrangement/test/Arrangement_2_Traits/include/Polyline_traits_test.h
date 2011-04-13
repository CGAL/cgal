#include <list>
#include "Base_traits_test.h"

template< class Traits_class, class Number_type >
class Polyline_traits_test : public Base_traits_test< Traits_class, Number_type >{
public:
  typedef Number_type  NT;
  typedef typename Traits_class::Point     Point;
  typedef typename Traits_class::X_curve   X_curve;
  typedef typename Traits_class::Curve     Curve;
public:
  Polyline_traits_test( int argc, char** argv );
  virtual void read_curve( std::ifstream& is, Curve& cv );
  virtual bool make_x_monotone_wrapper( std::istrstream& strLine );
  virtual bool curve_split_wrapper( std::istrstream& strLine );
  ~Polyline_traits_test();
};

//------------------------------------------------------------------------------
/*
  Constructor. Nothing to do. Just calls super class ctor
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
Polyline_traits_test< Traits_class, Number_type >::
Polyline_traits_test( int argc, char** argv ) :
 Base_traits_test< Traits_class, Number_type >(argc, argv) {}; 
//------------------------------------------------------------------------------
/*
  Destructor. Nothing to do. Implements super class virtual dtor.
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
Polyline_traits_test< Traits_class, Number_type >::
~Polyline_traits_test( ){}
//------------------------------------------------------------------------------
/*
  Reads on curve. This method is called by collect_data.
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
void Polyline_traits_test< Traits_class, Number_type >::
read_curve( std::ifstream& is, Curve& cv ){
  char one_line[128];
  int n_verteces;
  NT x,y;

  skip_comments( is, one_line );
  std::istrstream strLine( one_line, 128 );
  strLine >> n_verteces;
  for( int j = 0; j < n_verteces; j++ ){
    skip_comments( is, one_line );
    std::istrstream strLine( one_line, 128 );
    strLine >> x >> y;
    cv.push_back( Point( x,y ) );
  }
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
bool Polyline_traits_test< Traits_class, Number_type >::
make_x_monotone_wrapper( std::istrstream& strLine ){
  int index, exp_number, real_number;
  std::list<X_curve> l;
  typename std::list<X_curve>::iterator it;
  strLine >> index >> exp_number;
  
  tr.make_x_monotone( all_curves_vec[index], l );
  real_number = l.size();
  it = l.begin();
  Point prev_target_point = *( it->begin() );                   
  Point curr_source_point;
  std::cout << "Test: makes_x_monotone( Curve" << index << " ) ? " 
       << exp_number << std::endl;
  int nSubcurveIndex = 0;
  for( ; it != l.end(); ++it ){
    curr_source_point = *( it->begin() );
    if( !tr.is_x_monotone( *it ) ){
      std::cout << "Was NOT successful" << std::endl;
      std::cout << "One of the result subcurves is not x-monotone" << std::endl;
      std::cout << "Subcurve index: " << nSubcurveIndex << std::endl;
      return false;
    }
    if( curr_source_point != prev_target_point ){
      std::cout << "Was NOT successful" << std::endl;
      std::cout << "Some points are absent" << std::endl;
      std::cout << "Subcurve No " << nSubcurveIndex << "Source point: " 
                << curr_source_point << std::endl;
      std::cout << "Subcurve No " << nSubcurveIndex - 1 << "Target point: "
                << prev_target_point << std::endl;
      return false;
    }
    prev_target_point = *( it->rbegin() );                 // Assumption: container is bi-direcional
    ++nSubcurveIndex;
  }
  if( exp_number < real_number ){
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "Number of x-monotone subcurves is greater than expected " << std::endl;
    std::cout << "Expected: " << exp_number << "Actual: " << real_number << std::endl;
    return false;
  }
  else if( exp_number > real_number ){
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "Number of x-monotone subcurves is less than expected " << std::endl;
    std::cout << "Expected: " << exp_number << "Actual: " << real_number << std::endl;
    return false;  
  }
  else{
    std::cout << "Was successful" << std::endl;
    return true;
  }
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
bool Polyline_traits_test< Traits_class, Number_type >::
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
  if( *( (all_curves_vec[index1]).begin() )  == all_points_vec[index2] ||
      *( (all_curves_vec[index1]).rbegin() ) == all_points_vec[index2]  ){
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "Precondition fault: Split point is the end point of the curve " << std::endl;
    return false;
  }
  tr.curve_split( all_curves_vec[index1], cv1, cv2, all_points_vec[index2] );
  unsigned int i = 0,j = 0;
  for( ; i < cv1.size() - 1; i++ ){ //without the end point of cv1
    if( (all_curves_vec[index1])[ i ] != cv1[i] ){
      std::cout << "Was not successful " << std::endl;
      std::cout << "Some point is absent in the 1st part" << std::endl;
      std::cout << "Original curve[" << i << "]: " 
                << (all_curves_vec[index1])[ i ] << std::endl;
      std::cout << "obtained curve[ "<< i << "] " << cv1[i] << std::endl;          
      return false;
    }
  }
  if( ( all_curves_vec[index1])[i] != cv1[i] ){
    j = 1;
  }
  for( ; j < cv2.size(); j++ ){
    if( ( all_curves_vec[index1])[i] != cv2[j] ){
      std::cout << "Was not successful " << std::endl;
      std::cout << "Some point is absent in the 2nd part" << std::endl;
      std::cout << "Original curve[" << i << "]: " 
                << (all_curves_vec[index1])[ i ] << std::endl;
      std::cout << "obtained curve[ "<< j << "] " << cv2[j] << std::endl;          
      return false;
    }
    ++i;
  }
  std::cout << "Was successfull" << std::endl;
  return true;
}
//------------------------------------------------------------------------------
