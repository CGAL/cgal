// examples/Sweep_line/example4.C
// ------------------------------
#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h> 

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>

#include <CGAL/sweep_to_produce_planar_map_subcurves.h>

#include <CGAL/Arr_polyline_traits.h>
//#include <CGAL/IO/Arr_polyline_traits_iostream.h>

#include <iostream>
#include <vector>

typedef CGAL::Quotient<int>                  NT;
typedef CGAL::Cartesian<NT>                  R;
typedef CGAL::Arr_polyline_traits<R>         Traits;

typedef Traits::Point                                 Point;
typedef Traits::X_curve                               X_curve;
typedef Traits::Curve                                 Curve;

using namespace CGAL;


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
  std::cout<<"number of polylines is : " << num_polylines<<std::endl;

  while (num_polylines--) {
    Curve      polyline;
    
    std::cin>>polyline;
    
    curves.push_back(polyline);

    polyline.clear();
  }
}

int main()
{
  std::list<Curve>      polylines;
    
  read_polylines(polylines);
  
  Traits traits;
  std::list<Curve> subcurves;
  CGAL::sweep_to_produce_planar_map_subcurves(polylines.begin(),polylines.end(), traits, subcurves);  
 
  for (std::list<Curve>::iterator scv_iter = subcurves.begin(); scv_iter != subcurves.end(); scv_iter++)    
    std::cout<<*scv_iter<<std::endl;
}




