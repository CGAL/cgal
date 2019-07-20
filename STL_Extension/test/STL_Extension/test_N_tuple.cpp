#include <CGAL/internal/disable_deprecation_warnings_and_errors.h>

#include <CGAL/Twotuple.h>
#include <CGAL/Threetuple.h>
#include <CGAL/Fourtuple.h>
#include <CGAL/Sixtuple.h>
#include <CGAL/use.h>

int main()
{
  CGAL::Twotuple<int> d2, t2(0,1);
  CGAL::Threetuple<int> d3, t3(0,1,2);
  CGAL::Fourtuple<int> d4, t4(0,1,2,3);
  CGAL::Sixtuple<int> d6, t6(0,1,2,3,4,5);

  CGAL_USE(d2); CGAL_USE(t2);
  CGAL_USE(d3); CGAL_USE(t3);
  CGAL_USE(d4); CGAL_USE(t4);
  CGAL_USE(d6); CGAL_USE(t6);

  return 0;
}
