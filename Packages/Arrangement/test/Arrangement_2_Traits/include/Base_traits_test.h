#include <iostream>
#include <fstream>
#include <strstream>
#include <string>
#include <vector>

template< class T >
bool print_was_successful_or_not(  T& exp_answer, 
                                   T& real_answer );
void   skip_comments(              std::ifstream& is, 
                                   char* one_line );
std::string remove_blanks(         char* str );

template< class Traits_class, class Number_type >
class Base_traits_test{
public:
  typedef Number_type  NT;
  typedef typename Traits_class::Point     Point;
  typedef typename Traits_class::X_curve   X_curve;
  typedef typename Traits_class::Curve     Curve;

public:
  Base_traits_test( int argc, char** argv );
  bool start();
  virtual ~Base_traits_test();

protected:
  virtual void read_curve( std::ifstream& is, Curve& cv ) = 0;
  void collect_data( std::ifstream& is );
  bool perform_test( std::ifstream& is );
  bool compare_x_y_wrapper( std::istrstream& strLine, std::string& strCommand );
  bool curve_is_vertical_wrapper( std::istrstream& strLine );
  bool curve_is_in_x_range_wrapper( std::istrstream& strLine );
  bool curve_compare_at_x_smth_wrapper( std::istrstream& strLine, 
                                        std::string& strCommand  );
  bool curve_get_point_status_wrapper( std::istrstream& strLine );
  bool curve_is_between_cw_wrapper( std::istrstream& strLine );
  bool curve_is_same_wrapper( std::istrstream& strLine );
  bool curve_src_trg_wrapper( std::istrstream& strLine, std::string& strCommand );
  bool point_to_lr_wrapper( std::istrstream& strLine, std::string& strCommand );
  bool is_x_monotone_wrapper( std::istrstream& strLine );
  virtual bool make_x_monotone_wrapper( std::istrstream& strLine ) = 0;
  virtual bool curve_split_wrapper( std::istrstream& strLine ) = 0;
  bool point_reflect_in_x_and_y_wrapper( std::istrstream& strLine );
  bool curve_reflect_in_x_and_y_wrapper( std::istrstream& strLine );
  bool do_intersect_to_right_wrapper( std::istrstream& strLine );
  bool nearest_intersection_to_right_wrapper( std::istrstream& strLine );
  bool curves_overlap_wrapper( std::istrstream& strLine );

  bool get_expected_boolean( std::istrstream& strLine );
  int  get_expected_enum( std::istrstream& strLine );
  bool translate_boolean( std::string& strValue );
  int  translate_enumerator( std::string& strValue );

protected:
  std::string         data_file;
  Traits_class   tr;
  std::vector<Curve>  all_curves_vec;
  std::vector<Point>  all_points_vec;

};


//------------------------------------------------------------------------------
/*
  Constructor. 
  Recieves test data file name.
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
Base_traits_test< Traits_class, Number_type >::
Base_traits_test( int argc, char** argv ): tr(), 
                                           all_curves_vec(), 
                                           all_points_vec(){
  if( argc < 2 || argc >= 3 ){
    std::cout << "Usage: "<< argv[0] << " test_data_file" << std::endl;  
  }  
  else{
    data_file = argv[1];
  }
}
//------------------------------------------------------------------------------
/*
  Destructor. 
  Declares as virtual.
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
Base_traits_test< Traits_class, Number_type >::
~Base_traits_test( ){}
//------------------------------------------------------------------------------
/*
  Test entry point 
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
bool Base_traits_test< Traits_class, Number_type >::
start(){
  std::ifstream is( data_file.c_str() );
  if( !is ){
    return false;
  } 
  collect_data( is );
  bool test_result = perform_test( is );
  is.close(); 
  return test_result;
}
//------------------------------------------------------------------------------
/*
  Collects data. Fills all_curves_vec and all_points_vec
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
void Base_traits_test< Traits_class, Number_type >::
collect_data( std::ifstream& is ){
  char one_line[128];
  int  n_curves, n_points;
  NT x,y;
  int i;

  skip_comments( is, one_line );
  std::istrstream strLine( one_line, 128 );
  strLine >> n_curves;
  for( i = 0; i < n_curves; i++ ){
    Curve cv;
    read_curve( is, cv );
    all_curves_vec.push_back( cv );
  }
  skip_comments( is, one_line );
  std::istrstream strLine2( one_line, 128 );
  strLine2 >> n_points;
  for( i = 0; i < n_points; i++ ){
    skip_comments( is, one_line );
    std::istrstream strLine( one_line, 128 );
    strLine >> x >> y;
    all_points_vec.push_back( Point( x, y ) );
  } 
}
//------------------------------------------------------------------------------
/*
  Command dispatcher. Retrieves a line from the input file and performes
  some action. See comments for suitable function in order to know specific
  command arguments. 
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
bool Base_traits_test< Traits_class, Number_type >::
perform_test( std::ifstream& is ){
  char one_line[128];
  char buff[128];
  bool test_result = true; // initialize to goo test results
  std::cout << "Performing test ..." << std::endl;  
  while( !is.eof() ){
    buff[0] = '\0';
    skip_comments( is, one_line );
    std::istrstream strLine( one_line, 128 );
    strLine.getline( buff, 128, ' ' );
    std::string strCommand( buff );
    if( strCommand == "compare_x" || strCommand == "compare_y" ){
      test_result &= compare_x_y_wrapper( strLine, strCommand );
    }
    else if( strCommand == "curve_is_vertical" ){
      test_result &= curve_is_vertical_wrapper( strLine );
    }
    else if( strCommand == "curve_is_in_x_range" ){
      test_result &= curve_is_in_x_range_wrapper( strLine );
    }
    else if( strCommand == "curve_compare_at_x" ||
             strCommand == "curve_compare_at_x_left" ||
             strCommand == "curve_compare_at_x_right" ){
      test_result &= curve_compare_at_x_smth_wrapper( strLine, strCommand );
    }
    else if( strCommand == "curve_get_point_status" ){ 
      test_result &= curve_get_point_status_wrapper( strLine );
    }
    else if( strCommand == "curve_is_between_cw" ){
      test_result &= curve_is_between_cw_wrapper( strLine );
    }
    else if( strCommand == "curve_is_same" ){
      test_result &= curve_is_same_wrapper( strLine );
    }
    else if( strCommand == "curve_source" ||
             strCommand == "curve_target" ){
      test_result &= curve_src_trg_wrapper( strLine, strCommand );
    }
    else if( strCommand == "point_to_left"  ||
             strCommand == "point_to_right" ){
      test_result &= point_to_lr_wrapper( strLine, strCommand );
    }
    else if( strCommand == "is_x_monotone" ){
      test_result &= is_x_monotone_wrapper( strLine );
    }
    else if( strCommand == "make_x_monotone" ){
      test_result &= make_x_monotone_wrapper( strLine );
    }
    else if( strCommand == "curve_split" ){
      test_result &= curve_split_wrapper( strLine );
    }
    else if( strCommand == "point_reflect_in_x_and_y" ){
      test_result &= point_reflect_in_x_and_y_wrapper( strLine );
    }
    else if( strCommand == "curve_reflect_in_x_and_y" ){
      test_result &= curve_reflect_in_x_and_y_wrapper( strLine );
    }
    else if( strCommand == "do_intersect_to_right" ){
      test_result &= do_intersect_to_right_wrapper( strLine );
    }
    else if( strCommand == "nearest_intersection_to_right" ){
      test_result &= nearest_intersection_to_right_wrapper( strLine ); 
    }
    else if( strCommand == "curves_overlap" ){
      test_result &= curves_overlap_wrapper( strLine );
    }
  }

  return test_result;
}
//------------------------------------------------------------------------------
/*
  input case:
  compare_x n1, n2, RESULT or
  compare_y n1 n2 RESULT, where
  n1, n2 - points indeces in all_points_vec
  RESULT - expected result, enum{LARGER, SMALLER, EQUAL} type
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
bool Base_traits_test< Traits_class, Number_type >::
compare_x_y_wrapper( std::istrstream& strLine, std::string& strCommand ){
  int index1, index2, exp_answer, real_answer;

  strLine >> index1 >> index2;
  exp_answer = get_expected_enum( strLine );
  if( strCommand == "compare_x" ){
    real_answer = tr.compare_x( all_points_vec[index1], all_points_vec[index2] );
    std::cout << "Test: compare_x( " << all_points_vec[index1] 
         <<", "<< all_points_vec[index2] << " ) ? " << exp_answer << std::endl;
    return print_was_successful_or_not( exp_answer, real_answer );
  }
  else{
    real_answer = tr.compare_y( all_points_vec[index1], all_points_vec[index2] );
    std::cout << "Test: compare_y( " << all_points_vec[index1] 
         <<", "<< all_points_vec[index2] << " ) ? " << exp_answer << std::endl;
    return print_was_successful_or_not( exp_answer, real_answer );
  }
  return true;
}
//------------------------------------------------------------------------------
/*
  input case:
  curve_is_vertical n1 BOOL_RESULT, where
  n1 - curve index in all_curves_vec
  BOOL_RESULT - expected result, or FALSE or TRUE string 
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
bool Base_traits_test< Traits_class, Number_type >::
curve_is_vertical_wrapper( std::istrstream& strLine ){
  int index;
  bool exp_answer, real_answer;

  strLine >> index;
  exp_answer = get_expected_boolean( strLine );
  std::cout << "Test: curve_is_vertical( Curve" << index
       << " ) ? " << exp_answer << std::endl;
  if( !tr.is_x_monotone( all_curves_vec[index] ) ){
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "The input curve must be X-monotone" << std::endl;
    return false;    
  }
  real_answer = tr.curve_is_vertical( all_curves_vec[index] );
  return print_was_successful_or_not( exp_answer, real_answer );
}
//------------------------------------------------------------------------------
/*
  input case:
  curve_is_in_x_range n1 n2 BOOL_RESULT, where
  n1 - curve index in all_curves_vec
  n2 - point index in all_points_vec
  BOOL_RESULT - expected result, or FALSE or TRUE string
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
bool Base_traits_test< Traits_class, Number_type >::
curve_is_in_x_range_wrapper( std::istrstream& strLine ){
  int index1, index2;
  bool exp_answer, real_answer;

  strLine >> index1 >> index2;
  exp_answer = get_expected_boolean( strLine );
  std::cout << "Test: curve_is_in_x_range( Curve" << index1
       << ", " <<  all_points_vec[index2] << " ) ? " << exp_answer << std::endl;
  if( !tr.is_x_monotone( all_curves_vec[index1] ) ){
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "The input curve must be X-monotone" << std::endl;
    return false;    
  }
  real_answer = tr.curve_is_in_x_range( all_curves_vec[index1], 
                                        all_points_vec[index2] );
  return print_was_successful_or_not( exp_answer, real_answer );
}
//------------------------------------------------------------------------------
/*
  input case:
  curve_compare_at_x       n1 n2 n3 RESULT
  curve_compare_at_x_left  n1 n2 n3 RESULT
  curve_compare_at_x_right n1 n2 n3 RESULT, where
  n1 - curve index in all_curves_vec
  n2 - curve index in all_curves_vec
  n3 - point index in all_points_vec
  RESULT - expected result, enum{LARGER, SMALLER, EQUAL} type
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
bool Base_traits_test< Traits_class, Number_type >::
curve_compare_at_x_smth_wrapper( std::istrstream& strLine, 
                                 std::string& strCommand  ){
  int index1, index2, index3, exp_answer, real_answer;

  strLine >> index1 >> index2 >> index3;
  exp_answer = get_expected_enum( strLine );
  if( strCommand == "curve_compare_at_x_left" ){
    std::cout << "Test: curve_compare_at_x_left( Curve" << index1
         << ", Curve" << index2 << ", " 
         <<   all_points_vec[index3] << " ) ?" << std::endl 
	 << "  expected: " << exp_answer;
  }
  else if( strCommand == "curve_compare_at_x_right" ){
    std::cout << "Test: curve_compare_at_x_right( Curve" << index1
         << ", Curve" << index2 << ", " 
         <<   all_points_vec[index3] << " ) ?" << std::endl
	 << "  expected: " << exp_answer;
  }
  else{
    std::cout << "Test: curve_compare_at_x( Curve" << index1
         << ", Curve" << index2 << ", " 
         <<   all_points_vec[index3] << " ) ?" << std::endl
	 << "  expected: " << exp_answer;
  }
  if( ! ( tr.is_x_monotone( all_curves_vec[index1] )&&
          tr.is_x_monotone( all_curves_vec[index2] ) ) ){
    std::cout << std::endl << "  Was NOT successful" << std::endl;
    std::cout << "  The input curves must be X-monotone" << std::endl;
    return false;    
  }
  if( strCommand == "curve_compare_at_x_left" ){
    real_answer = tr.curve_compare_at_x_left( all_curves_vec[index1], 
                                              all_curves_vec[index2], 
                                              all_points_vec[index3] );

    std::cout << " actual: " << real_answer << std::endl;
  }
  else if( strCommand == "curve_compare_at_x_right" ){
    real_answer = tr.curve_compare_at_x_right( all_curves_vec[index1], 
                                               all_curves_vec[index2], 
                                               all_points_vec[index3] );
    std::cout << " actual: " << real_answer << std::endl;
  }
  else{
    real_answer = tr.curve_compare_at_x( all_curves_vec[index1], 
                                         all_curves_vec[index2], 
                                         all_points_vec[index3] );
    std::cout << " actual: " << real_answer << std::endl;
  }
  std::cout << "  ";
  return print_was_successful_or_not( exp_answer, real_answer );
}
//------------------------------------------------------------------------------
/*
  input case:
  curve_get_point_status n1 n2 STATUS_RESULT, where
  n1 - curve index in all_curves_vec
  n2 - point index in all_points_vec
  STATUS_RESULT - expected result, enum{ UNDER_CURVE, 
                                         CURVE_NOT_IN_RANGE, 
                                         ABOVE_CURVE, ON_CURVE } type
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
bool Base_traits_test< Traits_class, Number_type >::
curve_get_point_status_wrapper( std::istrstream& strLine ){
  int index1, index2, exp_answer, real_answer;

  strLine >> index1 >> index2;
  exp_answer = get_expected_enum( strLine );
  std::cout << "Test: curve_get_point_status( Curve" << index1
       << ", " <<  all_points_vec[index2] << " ) ? " << exp_answer << std::endl;
  if( !tr.is_x_monotone( all_curves_vec[index1] ) ){
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "The input curve must be X-monotone" << std::endl;
    return false;    
  }
  real_answer = tr.curve_get_point_status( all_curves_vec[index1], 
                                           all_points_vec[index2] );
  return print_was_successful_or_not( exp_answer, real_answer );
}
//------------------------------------------------------------------------------
/*
  input case:
  curve_is_between_cw n1 n2 n3 n4 BOOL_RESULT, where
  n1, n2, n3 - curves indeces in all_curves_vec
  n4 - point index in all_point_vec
  BOOL_RESULT - expected result, or FALSE or TRUE string  
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
bool Base_traits_test< Traits_class, Number_type >::
curve_is_between_cw_wrapper( std::istrstream& strLine ){
  int index1, index2, index3, index4;
  bool exp_answer, real_answer;

  strLine >> index1 >> index2 >> index3 >> index4;
  exp_answer = get_expected_boolean( strLine );
  std::cout << "Test: curve_is_between_cw( Curve" << index1
       << ", Curve"<< index2 << ", Curve" << index3 
       << ", " <<  all_points_vec[index4] << " ) ? " << exp_answer << std::endl;
  if( ! ( tr.is_x_monotone( all_curves_vec[index1] ) &&
          tr.is_x_monotone( all_curves_vec[index2] ) &&
          tr.is_x_monotone( all_curves_vec[index3] ) ) ){
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "The input curve must be X-monotone" << std::endl;
    return false;    
  }
  real_answer = tr.curve_is_between_cw( all_curves_vec[index1],
                                        all_curves_vec[index2],
                                        all_curves_vec[index3], 
                                        all_points_vec[index4] );
  return print_was_successful_or_not( exp_answer, real_answer );  
}
//------------------------------------------------------------------------------
/*
  input case:
  curve_is_same n1 n2  BOOL_RESULT, where
  n1, n2 - curves indeces in all_curves_vec
  BOOL_RESULT - expected result, or FALSE or TRUE string  
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
bool Base_traits_test< Traits_class, Number_type >::
curve_is_same_wrapper( std::istrstream& strLine ){
  int index1, index2, exp_answer, real_answer;

  strLine >> index1 >> index2;
  exp_answer = get_expected_boolean( strLine );
  std::cout << "Test: curve_is_same( Curve" << index1
       << ", Curve" <<  index2 << " ) ? " << exp_answer << std::endl;
  if( !( tr.is_x_monotone( all_curves_vec[index1] ) &&
         tr.is_x_monotone( all_curves_vec[index2] ) ) ){
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "The input curve must be X-monotone" << std::endl;
    return false;    
  }
  real_answer = tr.curve_is_same( all_curves_vec[index1], 
                                  all_curves_vec[index2] );
  return print_was_successful_or_not( exp_answer, real_answer );
}
//------------------------------------------------------------------------------
/*
  input case:
  curve_source  n1 x y
  curve_target  n1 x y, where
  n1 - curve index in all_curves_vec
  x, y - NT values, represent expected point X and Y coordinate
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
bool Base_traits_test< Traits_class, Number_type >::
curve_src_trg_wrapper( std::istrstream& strLine, std::string& strCommand ){
  int index1;
  Point real_answer;
  NT x,y;

  strLine >> index1 >> x >> y;
  Point exp_answer( x, y );
  if( strCommand == "curve_source" ){
    std::cout << "Test: curve_source( Curve" << index1 << ") ? " << exp_answer<<std::endl;
  }
  else{
    std::cout << "Test: curve_target( Curve" << index1 << ") ? " << exp_answer<<std::endl;
  }
  if( !tr.is_x_monotone( all_curves_vec[index1] ) ){
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "The input curve must be X-monotone" << std::endl;
    return false;    
  }
  if( strCommand == "curve_source" ){
    real_answer = tr.curve_source( all_curves_vec[index1] );
  }
  else{
    real_answer = tr.curve_target( all_curves_vec[index1] );
  }
  return print_was_successful_or_not( exp_answer, real_answer );
}
//------------------------------------------------------------------------------
/*
  input case:
  point_to_left  n1
  point_to_right n1
  n1 - point index in all_points_vec
  Does NOT take any expected result
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
bool Base_traits_test< Traits_class, Number_type >::
point_to_lr_wrapper( std::istrstream& strLine, std::string& strCommand ){
  int index;
  
  strLine >> index;
  if( strCommand == "point_to_left" ){
    std::cout << "Test: point_to_left( " << all_points_vec[index] << " )" << std::endl;
    Point answer = tr.point_to_left( all_points_vec[index] );
    if( tr.compare_x( answer, all_points_vec[index] ) == CGAL::SMALLER ){
      std::cout << "Was successful" << std::endl;
      return true;
    }
    else{
      std::cout << "Was NOT successful" << std::endl;
      std::cout << "Generated point is NOT to the left of the original." << std::endl;
      return false;
    }
  }
  else{
    std::cout << "Test: point_to_right( " << all_points_vec[index] << " )" << std::endl;
    Point answer = tr.point_to_right( all_points_vec[index] );
    if( tr.compare_x( answer, all_points_vec[index] ) == CGAL::LARGER ){
      std::cout << "Was successful" << std::endl;
      return true;
    }
    else{
      std::cout << "Was NOT successful" << std::endl;
      std::cout << "Generated point is NOT to the right of the original." << std::endl;
      return false;
    }
  }
}
//------------------------------------------------------------------------------
/*
  input case:
  is_x_monotone n1, BOOL_RESULT, where
  n1 - curve index in all_curves_vec
  BOOL_RESULT - expected result, or FALSE or TRUE string  
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
bool Base_traits_test< Traits_class, Number_type >::
is_x_monotone_wrapper( std::istrstream& strLine ){
  int index;
  bool exp_answer, real_answer;

  strLine >> index;
  exp_answer = get_expected_boolean( strLine );
  std::cout << "Test: is_x_monotone( Curve" << index
       << " ) ? " << exp_answer << std::endl;
  real_answer = tr.is_x_monotone( all_curves_vec[index] );
  return print_was_successful_or_not( exp_answer, real_answer );    
}
//------------------------------------------------------------------------------
/*
  input case:
  point_reflect_in_x_and_y  n1 x y, where
  n1 - point index in all_points_vec
  x, y - NT values, represent expected point X and Y coordinate
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
bool Base_traits_test< Traits_class, Number_type >::
point_reflect_in_x_and_y_wrapper( std::istrstream& strLine ){
  int index;
  NT exp_x, exp_y;

  strLine >> index >> exp_x >> exp_y;
  Point exp_answer( exp_x, exp_y );
  std::cout << "Test: point_reflect_in_x_and_y( " << all_points_vec[index] << " ) ?"
       << exp_answer << std::endl;
  Point real_answer = tr.point_reflect_in_x_and_y( all_points_vec[index] );
  return print_was_successful_or_not( exp_answer, real_answer );
}
//------------------------------------------------------------------------------
/*
  input case:
  curve_reflect_in_x_and_y n1, where 
  n1 - curve index in all_curve_vector
  Does NOT take any expected result
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
bool Base_traits_test< Traits_class, Number_type >::
curve_reflect_in_x_and_y_wrapper( std::istrstream& strLine ){
  int index;

  strLine >> index;
  std::cout << "Test: curve_reflect_x_and_y( Curve" << index << " )"
       << std::endl;
  if( ! tr.is_x_monotone( all_curves_vec[ index ] ) ){
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "Input curve is not x-monotone" << std::endl;
    return false;
  }
  X_curve res  = tr.curve_reflect_in_x_and_y( all_curves_vec[ index ] );
  X_curve res2 = tr.curve_reflect_in_x_and_y( res );
  if( tr.curve_is_same( res2, all_curves_vec[ index ] ) ){
    std::cout << "Was successful" << std::endl;
    return true;
  }
  else{
    std::cout << "Was NOT successful" << std::endl;
    return false;
  }
}
//------------------------------------------------------------------------------
/*
  input case:
  do_intersect_to_right n1 n2 n3 BOOL_RESULT, where
  n1, n2 - curves indeces in all_curves_vec
  n3 - point index in all_point_vec
  BOOL_RESULT - expected result, or FALSE or TRUE string  
  
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
bool Base_traits_test< Traits_class, Number_type >::
do_intersect_to_right_wrapper( std::istrstream& strLine ){
  int index1, index2, index3;
  bool exp_answer, real_answer;

  strLine >> index1 >> index2 >> index3;
  exp_answer = get_expected_boolean( strLine );
  std::cout << "Test: do_intersect_to_right( Curve" << index1 << ", "
       << "Curve " << index2 << ", " << all_points_vec[index3] 
       << " ) ? " << exp_answer << std::endl; 
  if( !( tr.is_x_monotone( all_curves_vec[index1] ) &&
         tr.is_x_monotone( all_curves_vec[index2] ) ) ){
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "The input curve must be X-monotone" << std::endl;
    return false;    
  }
  real_answer = tr.do_intersect_to_right( all_curves_vec[index1], 
                                          all_curves_vec[index2],
                                          all_points_vec[index3] );
  return print_was_successful_or_not( exp_answer, real_answer );
}
//------------------------------------------------------------------------------
/*
  input case:
  nearest_intersection_to_right n1 n2 n3 x1 y1 x2 y2, where
  n1, n2 - curves indeces in all_curves_vec
  n3 - point index in all_point_vec
  x1, y1 - NT values, represent expected intersection 
                                                SOURCE point X and Y coordinate
  x2, y2 - NT values, represent expected intersection 
                                               TARGET point X and Y coordinate
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
bool Base_traits_test< Traits_class, Number_type >::
nearest_intersection_to_right_wrapper( std::istrstream& strLine ){
  int index1, index2, index3;
  NT x1, y1, x2, y2;
  bool exp_answer, real_answer;
  Point p1, p2;

  strLine >> index1 >> index2 >> index3 >> x1 >> y1 >> x2 >> y2;
  Point exp_p1( x1, y1 ), exp_p2( x2, y2 );
  exp_answer = get_expected_boolean( strLine );
  std::cout << "Test: nearest_intersection_to_right( Curve" << index1 << ", "
       << "Curve " << index2 << ", " << all_points_vec[index3] 
       << " ) ? " << exp_answer << std::endl; 
  if( !( tr.is_x_monotone( all_curves_vec[index1] ) &&
         tr.is_x_monotone( all_curves_vec[index2] ) ) ){
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "The input curve must be X-monotone" << std::endl;
    return false;    
  }
  real_answer = tr.nearest_intersection_to_right( all_curves_vec[index1],
                                                  all_curves_vec[index2],
                                                  all_points_vec[index3],
                                                  p1, p2 );
  if( exp_answer == real_answer && p1 == exp_p1 && p2 == exp_p2 ){
    std::cout << "Was successful" << std::endl;
    return true;
  }
  else if( exp_answer == real_answer && exp_answer == false ){
    // We were expecting to get FALSE.
    // So it does not matter what are the points.
    std::cout << "Was successful" << std::endl;
    return true;
  }
  else{
    std::cout << "Was NOT successful" << std::endl;
    if( exp_answer != real_answer ){
      std::cout << "The returned answer is different from the expected one" << std::endl;
    }
    if( p1 != exp_p1 ){
      std::cout << "Source point of the intersection is not the expected one" << std::endl;
      std::cout << "Expected: " << exp_p1 << " Actual: " << p1 << std::endl;
    }
    if( p2 != exp_p2 ){
      std::cout << "Target point of the intersection is not the expected one" << std::endl;
      std::cout << "Expected: " << exp_p2 << " Actual: " << p2 << std::endl;
    }
    return false;
  }
} 
//------------------------------------------------------------------------------
/*
  input case:
  curves_overlap n1 n2 BOOL_RESULT, where
  n1, n2 - curves indeces in all_curves_vec
  BOOL_RESULT - expected result, or FALSE or TRUE string  
*/
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
bool Base_traits_test< Traits_class, Number_type >::
curves_overlap_wrapper( std::istrstream& strLine ){
  int index1, index2;
  bool exp_answer, real_answer;

  strLine >> index1 >> index2;
  exp_answer = get_expected_boolean( strLine );
  std::cout << "Test: curves_overlap( Curve" << index1 << ", "
       << "Curve " << index2
       << " ) ? " << exp_answer << std::endl; 
  if( !( tr.is_x_monotone( all_curves_vec[index1] ) &&
         tr.is_x_monotone( all_curves_vec[index2] ) ) ){
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "The input curve must be X-monotone" << std::endl;
    return false;    
  }
  real_answer = tr.curves_overlap( all_curves_vec[index1],
                                   all_curves_vec[index2] );
  return print_was_successful_or_not( exp_answer, real_answer );    
  
}
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
bool Base_traits_test< Traits_class, Number_type >::
translate_boolean( std::string& strValue ){
  if(  strValue == "TRUE" ){
    return true;
  }
  return false;
}
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
int Base_traits_test< Traits_class, Number_type >::
translate_enumerator( std::string& strValue ){
  if( strValue == "LARGER" ){
    return CGAL::LARGER;
  }
  else if( strValue == "SMALLER" ){
    return CGAL::SMALLER;
  }
  else if( strValue == "EQUAL" ){
    return CGAL::EQUAL;
  }
  else if( strValue == "UNDER_CURVE" ){ 
    return tr.UNDER_CURVE;
  }
  else if( strValue == "CURVE_NOT_IN_RANGE" ){
    return tr.CURVE_NOT_IN_RANGE;
  }
  else if( strValue == "ABOVE_CURVE" ){
    return tr.ABOVE_CURVE;
  }
  else if( strValue == "ON_CURVE" ){
    return tr.ON_CURVE;
  }

  return -220776; // My birthday :-)
}
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
bool Base_traits_test< Traits_class, Number_type >::
get_expected_boolean( std::istrstream& strLine ){
  char buff[128];
  strLine.getline( buff, 128, '.' );
  buff[strLine.gcount()] = '\0';
  std::string strExpRes = remove_blanks( buff );
  return translate_boolean( strExpRes );
}
//------------------------------------------------------------------------------
template< class Traits_class, class Number_type >
int Base_traits_test< Traits_class, Number_type >::
get_expected_enum( std::istrstream& strLine ){
  char buff[128];
  strLine.getline( buff, 128, '.' );
  buff[strLine.gcount()] = '\0';
  std::string strExpRes = remove_blanks( buff );
  return translate_enumerator( strExpRes );
}

//------------------------------------------------------------------------------
std::string remove_blanks( char* str ){
  std::string result = "";
  while( *str != '\0' ){
    if( *str != ' ' )
      result += *str;
    ++str;
  }
  return result;
}
//------------------------------------------------------------------------------
template< class T >
bool print_was_successful_or_not( T& exp_answer, T& real_answer ){
  bool bReturnValue = true;
  if( exp_answer == real_answer ){
    std::cout << "Was successful"     << std::endl;
  }
  else{
    std::cout << "Was NOT successful" << std::endl;
    std::cout << "Expected result: " << exp_answer << std::endl;
    std::cout << "Obtained result: " << real_answer << std::endl;
    bReturnValue = false;
  } 
  return bReturnValue;      
}
//------------------------------------------------------------------------------
/*
  Keeps on reading file lines while they are comments, i.e. starts with #-symbol
*/
//------------------------------------------------------------------------------
void skip_comments( std::ifstream& is, char* one_line ){
  while( !is.eof() ){
    is.getline( one_line, 128 );
    if( one_line[0] != '#' ){
      break;
    }
  }  
}
//------------------------------------------------------------------------------


