//examples/Arrangement_2/example12

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h> 

#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arrangement_2.h>

#include <iostream>
#include <CGAL/Arr_polyline_traits.h>
#include <CGAL/IO/Arr_iostream.h>


typedef CGAL::Quotient<int>                            NT;
typedef CGAL::Cartesian<NT>                            R;
typedef CGAL::Arr_polyline_traits<R>                   Traits;

typedef Traits::Point                                  Point;
typedef Traits::X_curve                                X_curve;
typedef Traits::Curve                                  Curve;

typedef CGAL::Arr_2_default_dcel<Traits>               Dcel;
typedef CGAL::Arr_base_node<Curve>                     Base_node;
typedef CGAL::Arrangement_2<Dcel,Traits,Base_node >    Arrangement;

typedef CGAL::Arr_file_writer<Arrangement>             Arr_writer;


CGAL_BEGIN_NAMESPACE

std::ostream&  operator<<(std::ostream& os,  
			  const Curve& cv)
{
  typedef Curve::const_iterator Points_iterator;
  
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
  unsigned int i;

  in >> size;

  for (i = 0; i < size; i++){
    Point  p;
    
    in >> p;
    
    cv.push_back(p);  
  }
  
  return in;
}

CGAL_END_NAMESPACE

int main()
{
  Arrangement arr;
  
  std::cin>>arr;
 
  std::cout << " * * * Printing list of all halfedges " 
	    << "of the resulting Arrangement" 
	    << std::endl;
  
  Arr_writer verbose_writer(std::cout, arr, true);
  
  verbose_writer.write_halfedges(arr.halfedges_begin(), arr.halfedges_end());

  return 0;
}
