/*
test the different convex_hull algorithm implementations
with the kernel archetype ...
*/

#include <CGAL/basic.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/ch_akl_toussaint.h>
#include <CGAL/ch_bykat.h>
#include <CGAL/ch_eddy.h>
#include <CGAL/ch_graham_andrew.h>
#include <CGAL/ch_jarvis.h>
#include <CGAL/ch_melkman.h>
#include <CGAL/Kernel_archetype.h>
#include <vector>
#include <list>

typedef CGAL::Kernel_archetype  K;
typedef K::Point_2              Point_2;

int main()
{
  std::vector<Point_2> input;
  
  Point_2 act;
  input.push_back(act);

  std::list<Point_2> output;

  K  traits;

  CGAL::convex_hull_2(input.begin(), input.end(), std::back_inserter(output), traits);
  
  CGAL::ch_akl_toussaint(input.begin(), input.end(), std::back_inserter(output), traits);
  
  CGAL::ch_bykat(input.begin(), input.end(), std::back_inserter(output), traits);
  
  CGAL::ch_eddy(input.begin(), input.end(), std::back_inserter(output), traits);
  
  CGAL::ch_graham_andrew(input.begin(), input.end(), std::back_inserter(output), traits);
  
  CGAL::ch_jarvis(input.begin(), input.end(), std::back_inserter(output), traits);
  
  CGAL::ch_melkman(input.begin(), input.end(), std::back_inserter(output), traits);            		       
  
  return 0;
}
