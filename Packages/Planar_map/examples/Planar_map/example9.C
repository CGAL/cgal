// example9

#include <CGAL/Cartesian.h>

//#include <CGAL/leda_rational.h>
#include <CGAL/Quotient.h>

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_segment_exact_traits.h>

#include <CGAL/IO/Pm_iostream.h>
#include <CGAL/IO/write_pm.h>
#include <iostream.h>

//uncomment if you have LEDA installed.
//#include <CGAL/IO/Pm_Window_stream.h>
//#include <CGAL/IO/leda_window.h>

//uncomment if you have LEDA installed.
//#include <CGAL/IO/Pm_Postscript_file_stream.h>

//using std::cout;
//using std::cin;
//using std::endl;
//using CGAL::write_pm;

//typedef leda_rational                    NT;
typedef CGAL::Quotient<int>                NT;
typedef CGAL::Cartesian<NT>                R;
typedef CGAL::Pm_segment_exact_traits<R>   Traits;
typedef CGAL::Pm_default_dcel<Traits>      Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>    PM;
typedef CGAL::Pm_file_writer<PM>           Pm_writer;

int main()
{ 
  PM pm;

  std::cout << "* * * Demonstrating a trivial use of IO functions" << std::endl << std::endl;
  std::cin  >> pm;
  std::cout << pm;
  
  std::cout << std::endl;
  std::cout << "* * * Presenting the use of verbose format" << std::endl;
  std::cout << std::endl;
  Pm_writer verbose_writer(std::cout, pm, true);
  write_pm(pm, verbose_writer, std::cout);
  
  std::cout << std::endl;
  std::cout << "* * * Demonstrating the use of the writer class interface." << std::endl;
  std::cout << "* * * Printing all halfedges in non verbose format" << std::endl << std::endl;
  Pm_writer writer(std::cout, pm);
  writer.write_halfedges(pm.halfedges_begin(), pm.halfedges_end());
  std::cout << std::endl;
  std::cout << "* * * Printing all halfedges in a verbose format" << std::endl << std::endl;
  verbose_writer.write_halfedges(pm.halfedges_begin(), pm.halfedges_end());
   
  //uncomment if you have LEDA installed. 
  //printing to leda window.
  //CGAL::Window_stream W(800, 800);
  //W.init(-5, +5, -5);
  //W.set_mode(leda_src_mode);
  //W.set_node_width(3);
  //W.display();
  //W << pm;
  
  //uncomment if you have LEDA installed.
  // printing to Postscript file.
  //CGAL::Postscript_file_stream  LPF(500, 500 ,"pm.ps");
  //LPF.init(-3,3,-3);
  //LPF.set_line_width( 1);
  //LPF << pm;
  
  return 0;
}
