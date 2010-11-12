#include <CGAL/Simple_cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Cartesian_converter.h>

typedef CGAL::Simple_cartesian<double>                           K1;
typedef CGAL::Simple_cartesian<CGAL::Quotient<CGAL::MP_Float> >  K2;
typedef CGAL::Cartesian_converter<K1,K2>                         K1_2_K2;
typedef CGAL::Cartesian_converter<K2,K1>                         K2_2_K1;

int main(){
  K1::Triangle_3 t_i(
    K1::Point_3(0.,0.,0.),
    K1::Point_3(1.,0.,-1.),
    K1::Point_3(0.,1.,3.)
  );
  
  K1::Line_3 l_i(
    K1::Point_3(0.2,0.25,-7),
    K1::Point_3(0.25,0.3,4)
  );
  
  K1_2_K2 to_exact;
  
  K2::Triangle_3 t_1=to_exact(t_i);
  K2::Line_3     l_e=to_exact(l_i);
  
  CGAL::Object inter=CGAL::intersection(t_e,l_e);
  K2::Point_3 exact_pt=CGAL::object_cast<K2::Point_3>(inter);
  
  K2_2_K1 to_inexact;
  
  K1::Point_3 approx_pt = to_inexact(exact_pt);
  return 0;
}
