// example9

#include <CGAL/basic.h> //CGAL definitions that need to come before anything

//#include <CGAL/leda_rational.h>
#include <CGAL/Quotient.h>

#include <CGAL/Cartesian.h>

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_segment_exact_traits.h>

#include <CGAL/IO/Pm_iostream.h>
#include <iostream.h>

//uncomment if you have LEDA installed.
//#include <CGAL/IO/Pm_Window_stream.h>
//#include <CGAL/IO/leda_window.h>

//uncomment if you have LEDA installed.
//#include <CGAL/IO/Pm_Postscript_file_stream.h>

//typedef leda_rational                    NT;
typedef CGAL::Quotient<int>                NT;
typedef CGAL::Cartesian<NT>                R;
typedef CGAL::Pm_segment_exact_traits<R>   Traits;
typedef CGAL::Pm_default_dcel<Traits>      Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>    PM;

int main()
{ 
  PM pm;

  std::cout<<"Demonstrating a trivial use of IO functions"<<std::endl;
  std::cin >> pm;
  std::cout << pm;
    
  std::cout<<"Presenting the use of verbose format"<<std::endl;
  CGAL::Pm_file_writer<PM>  verbose_writer(cout, pm, true);
  write_pm(pm, verbose_writer, cout);
  
  std::cout<<"Demonstrating the use of the writer class interface."<<std::endl;
  std::cout<<"Printing all halfedges in non verbose format"<<std::endl;
  CGAL::Pm_file_writer<PM>  writer(cout, pm);
  writer.write_halfedges(pm.halfedges_begin(), pm.halfedges_end());
  std::cout<<"Printing all halfedges in a verbose format"<<std::endl;
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

