#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Float_snap_rounding_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel         Kernel;
typedef Kernel::Segment_2                        Segment_2;
typedef Kernel::Point_2                          Point_2;
typedef Kernel::FT                               FT;

int main()
{
  std::vector<Point_2> pts;
  std::vector< std::vector<size_t> > segs;

  Kernel::FT eps(std::pow(2, -60));
  Kernel::FT one(1.);
  Kernel::FT two(2.);

  pts.emplace_back(one-eps, one);
  pts.emplace_back(-one+eps, -one+2*eps);
  pts.emplace_back(eps/2, eps/2);
  pts.emplace_back(one, -one);
  pts.emplace_back(0, two-eps/2);
  pts.emplace_back(two, 0);
  pts.emplace_back(-two+eps, -FT(4.));
  pts.emplace_back(two, two);
  pts.emplace_back(-two, two);

  segs.push_back(std::vector<size_t>());
  segs[0].push_back(0);
  segs[0].push_back(1);
  segs.push_back(std::vector<size_t>());
  segs[1].push_back(2);
  segs[1].push_back(3);
  segs.push_back(std::vector<size_t>());
  segs[2].push_back(4);
  segs[2].push_back(5);
  segs.push_back(std::vector<size_t>());
  segs[3].push_back(4);
  segs[3].push_back(6);
  segs.push_back(std::vector<size_t>());
  segs[4].push_back(7);
  segs[4].push_back(8);

  std::cout << "Input" << std::endl;
  std::cout << "Points" << std::endl;
  for(auto &p: pts)
    std::cout << p << std::endl;
  std::cout << "Segments:" << std::endl;
  for(auto &seg: segs){
    for(auto pi : seg)
        std::cout << pi << " ";
    std::cout << std::endl;
  }
  std::cout << "\n\n" << std::endl;

  CGAL::float_snap_rounding_2(pts, segs);

  std::cout << "Output" << std::endl;
  std::cout << "Points" << std::endl;
  for(auto &p: pts)
    std::cout << p << std::endl;
  std::cout << "Segments:" << std::endl;
  for(auto &seg: segs){
    for(auto pi : seg)
        std::cout << pi << " ";
    std::cout << std::endl;
  }

  return(0);
}
