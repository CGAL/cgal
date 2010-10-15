#include <CGAL/Complexity_tags.h>
#include <CGAL/Location_policy.h>

int main()
{
  CGAL::Fast f;
  CGAL::Compact c;
  CGAL::Fast_location fl;
  CGAL::Compact_location cl;
  CGAL::Location_policy<CGAL::Fast> flp = fl;
  CGAL::Location_policy<CGAL::Compact> clp = cl;

  (void) f;
  (void) c;
  (void) fl;
  (void) cl;
  (void) flp;
  (void) clp;
}
