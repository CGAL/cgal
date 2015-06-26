
// Without TBB_USE_THREADING_TOOL Intel Inspector XE will report false positives in Intel TBB
// (http://software.intel.com/en-us/articles/compiler-settings-for-threading-error-analysis-in-intel-inspector-xe/)
#ifdef _DEBUG
# define TBB_USE_THREADING_TOOL
#endif

#include <CGAL/Epick_d.h>
#include <CGAL/Tangential_complex/utilities.h>

#include <vector>

using namespace CGAL::Tangential_complex_;

void test_does_voronoi_face_and_alpha_tangent_subspace_intersect()
{
  typedef CGAL::Epick_d<CGAL::Dimension_tag<2> > K;
  typedef K::Point_d Pt;
  typedef K::Vector_d Vec;
  std::vector<std::size_t> P;
  std::vector<std::size_t> Q;
  std::vector<Vec> osb;
    
  K k;

  std::vector<Pt> points;
  points.push_back(Pt(0.02, -0.03));
  points.push_back(Pt(0.005, 2.3));
  points.push_back(Pt(4.5, 1.12));
  points.push_back(Pt(-3.5, 1.02));

  P.push_back(0);
  P.push_back(1);
    
  Q.push_back(2);
  Q.push_back(3);

  osb.push_back(Vec(0.01, 0.995));

  assert(does_voronoi_face_and_fixed_alpha_tangent_subspace_intersect(
    points, 0, P, Q, osb, 0.0, k) == false);
    
  assert(does_voronoi_face_and_fixed_alpha_tangent_subspace_intersect(
    points, 0, P, Q, osb, 0.5, k) == false);
    
  assert(does_voronoi_face_and_fixed_alpha_tangent_subspace_intersect(
    points, 0, P, Q, osb, 1.0, k) == false);
    
  assert(does_voronoi_face_and_fixed_alpha_tangent_subspace_intersect(
    points, 0, P, Q, osb, 1.5, k) == true);
}

int main()
{
  test_does_voronoi_face_and_alpha_tangent_subspace_intersect();
  return 0;
}
