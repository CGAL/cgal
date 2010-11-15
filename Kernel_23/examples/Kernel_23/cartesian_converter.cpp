#include <CGAL/Simple_cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Cartesian_converter.h>

typedef CGAL::Simple_cartesian<double>                           IK;
typedef CGAL::Simple_cartesian<CGAL::Quotient<CGAL::MP_Float> >  EK;
typedef CGAL::Cartesian_converter<IK,EK>                         IK_to_EK;
typedef CGAL::Cartesian_converter<EK,IK>                         EK_to_IK;

int main(){
  IK::Triangle_3 t1(
    IK::Point_3(0.,0.,0.),
    IK::Point_3(1.,0.,-1.),
    IK::Point_3(0.,1.,3.)
  );
  
  IK::Line_3 l1(
    IK::Point_3(0.2,0.25,-7),
    IK::Point_3(0.25,0.3,4)
  );
  
  IK_to_EK to_exact;
  
  EK::Triangle_3 t2=to_exact(t1);
  EK::Line_3     l2=to_exact(l1);
  
  CGAL::Object inter=CGAL::intersection(t2,l2);
  const EK::Point_3& exact_pt=CGAL::object_cast<EK::Point_3>(inter);
  
  EK_to_IK to_inexact;
  
  IK::Point_3 inexact_pt = to_inexact(exact_pt);
  return 0;
}
