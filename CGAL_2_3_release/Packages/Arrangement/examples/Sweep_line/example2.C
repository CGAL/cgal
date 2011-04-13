// examples/Sweep_line/example2.C
// ------------------------------
#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h> 

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>

#include <CGAL/sweep_to_construct_planar_map.h>

#include <CGAL/Arr_polyline_traits.h>

#include <CGAL/IO/Pm_iostream.h>

#include <iostream>
#include <vector>

//uncomment if you have LEDA installed.
//#ifndef CGAL_IO_PM_WINDOW_STREAM_H
//#include <CGAL/IO/Pm_Window_stream.h>
//#endif
//#ifndef CGAL_ARR_POLYLINE_TRAITS_WINDOW_STREAM_H  
//#include <CGAL/IO/Arr_polyline_traits_Window_stream.h>
//#endif
//#include <CGAL/IO/leda_window.h>  //used for visualization.

using namespace CGAL;


typedef CGAL::Quotient<int>                  NT;
typedef CGAL::Cartesian<NT>                  R;
typedef CGAL::Arr_polyline_traits<R>         Traits;

typedef Traits::Point                                 Point;
typedef Traits::X_curve                               X_curve;
typedef Traits::Curve                                 Curve;

typedef CGAL::Pm_default_dcel<Traits>                Dcel;   
typedef CGAL::Planar_map_2<Dcel, Traits>             PM;
typedef CGAL::Pm_file_writer<PM>                     Pm_writer;


CGAL_BEGIN_NAMESPACE

std::ostream&  operator<<(std::ostream& os,  
			  const Curve& cv)
{
  typedef Curve::const_iterator       Points_iterator;
  
  os<<cv.size()<<std::endl;
  for (Points_iterator points_iter = cv.begin(); 
       points_iter != cv.end(); points_iter++)
    os<<" "<<*points_iter;

  return os;
}


std::istream&  operator>>(std::istream& in,  
			  Curve& cv)
{
  std::size_t  size;

  in >> size;

  for (unsigned int i = 0; i < size; i++){
    Point  p;
    
    in >> p;
    
    cv.push_back(p);  
  }
  
  return in;
}

CGAL_END_NAMESPACE

template <class Container>
void read_polylines(Container& curves)
{
  int      num_polylines = 0;

  std::cin >> num_polylines;
  std::cout << "number of polylines is : " << num_polylines << std::endl;

  while (num_polylines--) {
    Curve      polyline;
    
    std::cin>>polyline;
    curves.push_back(polyline);
  }
}

int main()
{
  PM                 pm;
  std::vector<Curve>      polylines;
  
  read_polylines(polylines);

  CGAL::sweep_to_construct_planar_map(polylines.begin(),polylines.end(), pm);
  
  std::cout << " * * * Printing list of all halfedges of the resulting Planar map" << std::endl;
  
  Pm_writer verbose_writer(std::cout, pm, true);
  
  verbose_writer.write_halfedges(pm.halfedges_begin(), pm.halfedges_end());

  
  //CGAL::Window_stream W(700, 700);
  //W.init(-10, 10, -10);
  //W.set_mode(leda_src_mode);
  //W.set_node_width(3);
  //W.button("finish",2);
  //W.display();
  //W << pm;
}





