#include <CGAL/basic.h>
#include <CEP/Vanilla/Flavored_object.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Circle_2.h>

typedef double                     NT;
typedef CGAL::Cartesian<NT>        R;
typedef CGAL::Circle_2<R>          Circle;
typedef CGAL::Point_2<R>           Center;
typedef Flavored_object<Circle>    Scoop;

int main(int argc, char** argv)
{
   Scoop   ice_cream;
   ice_cream.set_flavor(PISTACHIO);

   if (! ice_cream.is_valid() || ice_cream.flavor() != PISTACHIO)
      exit (1);

   ice_cream.enhance_flavor();
   if (! ice_cream.is_valid() || ice_cream.flavor() != VANILLA )
      exit(1);

   Scoop   vanilla_scoop = Scoop(Circle(Center(3,2),8));

   if (!vanilla_scoop.is_valid() || vanilla_scoop.flavor() != VANILLA)
      exit (1);
   
   Scoop   chocolate_scoop(CHOCOLATE);
   if (!chocolate_scoop.is_valid() || chocolate_scoop.flavor() != CHOCOLATE)
      exit (1);

   Scoop   peach = Scoop(Circle(Center(0,0),16),PEACH);
   if (!peach.is_valid() || peach.flavor() != PEACH)
      exit (1);

   exit (0);
}
