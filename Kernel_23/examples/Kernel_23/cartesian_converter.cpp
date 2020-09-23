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

  CGAL::cpp11::result_of<EK::Intersect_3(EK::Triangle_3, EK::Line_3)>::type
    inter = CGAL::intersection(t2,l2);

  // As we are sure that there IS an intersection
  // and that the intersection IS a point
  // we do not have to check for this, or put it in a try/catch
  const EK::Point_3& exact_pt = boost::get<EK::Point_3>(*inter);

  EK_to_IK to_inexact;

  IK::Point_3 inexact_pt = to_inexact(exact_pt);
  std::cout << inexact_pt << std::endl;
  return 0;
}
