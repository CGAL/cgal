// example11

#include <CGAL/basic.h> //CGAL definitions that need to come before anything

#include <CGAL/Quotient.h>

#include <CGAL/Cartesian.h>

#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arr_segment_exact_traits.h>
#include <CGAL/Arrangement_2.h>

#include <CGAL/IO/Arr_iostream.h>
#include <iostream.h>

//uncomment if you have LEDA installed.
//#include <CGAL/IO/Arr_Window_stream.h>
//#include <CGAL/IO/leda_window.h>

//uncomment if you have LEDA installed.
//#include <CGAL/IO/Arr_Postscript_file_stream.h>

typedef CGAL::Quotient<int>                            NT;
typedef CGAL::Cartesian<NT>                            R;
typedef CGAL::Arr_segment_exact_traits<R>              Traits;

typedef Traits::Point                                 Point;
typedef Traits::X_curve                               X_curve;
typedef Traits::Curve                                 Curve;

typedef CGAL::Arr_2_default_dcel<Traits>               Dcel;
typedef CGAL::Arr_base_node<Curve>                     Base_node;
typedef CGAL::Arrangement_2<Dcel,Traits,Base_node >    Arrangement;

int main()
{ 
  Arrangement arr;

  std::cout<<"* * * Demonstrating a trivial use of IO functions"<<std::endl<<std::endl;
  std::cin >> arr;
  std::cout << arr;
  
  std::cout<<std::endl;

  std::cout<<"* * * Presenting the use of verbose format"<<std::endl<<std::endl;;
  CGAL::Arr_file_writer<Arrangement>  verbose_writer(cout, arr, true);
  write_arr(arr, verbose_writer, cout);
   
  //uncomment if you have LEDA installed. 
  //printing to leda window.
  //CGAL::Window_stream W(800, 800);
  //W.init(-5, +5, -5);
  //W.set_mode(leda_src_mode);
  //W.set_node_width(3);
  //W.display();
  //W << arr;
  
  //uncomment if you have LEDA installed.
  // printing to Postscript file.
  //CGAL::Postscript_file_stream  LPF(500, 500 ,"pm.ps");
  //LPF.init(-3,3,-3);
  //LPF.set_line_width( 1);
  //LPF << arr;
  
  return 0;
}
