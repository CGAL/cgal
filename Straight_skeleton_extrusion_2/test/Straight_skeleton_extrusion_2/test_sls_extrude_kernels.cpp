#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Exact_rational.h>
#ifdef CGAL_USE_GMPXX
#include <CGAL/gmpxx.h>
#endif
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>

#include "sls_extrusion_test_helpers.h"

#include <CGAL/extrude_skeleton.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <iostream>
#include <string>
#include <vector>

using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;
using EPECK_with_sqrt = CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt;
using SC_rational = CGAL::Simple_cartesian<CGAL::Exact_rational>;

namespace SSE = CGAL::Straight_skeleton_extrusion;
namespace SST = CGAL::Straight_skeletons_2::Extrusion_tests;

struct Case_description
{
  std::string name;
  double volume; // closed form
  double exact_num; // for exact kernels, volume = exact_num / exact_den
  double exact_den;
};

template <typename K>
std::vector<std::pair<Case_description, SST::Extrusion_input<K> > > kernel_cases()
{
  using FT = typename K::FT;
  using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

  using Input = SST::Extrusion_input<K>;
  using Weights = typename Input::Weights;

  std::vector<std::pair<Case_description, Input> > cases;

  cases.push_back({{"square pyramid", 4. / 3., 4, 3},
                  Input(SST::make_rectangle<K>(0, 0, 2, 2))});
  cases.push_back({{"rectangle hip roof", 10. / 3., 10, 3},
                  Input(SST::make_rectangle<K>(0, 0, 4, 2))});
  cases.push_back({{"L-shape", 16. / 3., 16, 3},
                  Input(SST::make_polygon<K>({{0,0}, {4,0}, {4,2}, {2,2}, {2,4}, {0,4}}))});
  cases.push_back({{"plus", 28. / 3., 28, 3},
                  Input(SST::make_polygon<K>({{2,0}, {4,0}, {4,2}, {6,2}, {6,4}, {4,4},
                                              {4,6}, {2,6}, {2,4}, {0,4}, {0,2}, {2,2}}))});
  {
    Input in(SST::make_rectangle<K>(0, 0, 2, 2));
    in.with_height(FT(1) / FT(2));
    cases.push_back({{"square frustum H = 1/2", 7. / 6., 7, 6}, in});  // 4/3 - 1/6
  }

  {
    Polygon_with_holes_2 ring(SST::make_rectangle<K>(0, 0, 10, 10));
    ring.add_hole(SST::make_rectangle<K>(4, 4, 6, 6, false));
    // decomposed into four pyramids + four "slices" + 8 tetrahedra
    // (4*4*2)*4/3 + (2*4*4) + 8*8/3
    cases.push_back({{"square ring", 96., 96, 1}, Input(ring)});

    Input in(ring);
    in.with_height(1);
    cases.push_back({{"square ring, H = 1", 72., 72, 1}, in}); // 96 - (2*2*1/3)*4 + (2*4*2/4)*4 + 8 * (1/3)

    Input out(ring);
    out.with_weights(Input::uniform(ring, -1.)).with_height(FT(1) / FT(2));
    cases.push_back({{"outward ring, H = 1/2", 54., 54, 1}, out}); // (5*5*2.5/3)*4 + (1*5*5/4)*4 + 8 * (12.5 / 2 * 2.5 / 3) - 96
  }

  {
    Input in(SST::make_rectangle<K>(0, 0, 4, 2));
    in.with_weights(Weights{{FT(0.5), FT(1), FT(0.5), FT(1)}});
    cases.push_back({{"weighted simultaneous event", 8. / 3., 8, 3}, in});
  }

  {
    // non-axis-aligned edges: constructions involve square roots
    Input in(SST::make_polygon<K>({{0,3}, {-1,1}, {-3,1}, {-1.5,-0.5}, {-2,-3}, {0,-1.5}, {2,-3}, {1.5,-0.5}, {3,1}, {1,1}}));
    // no idea what the volume is (but this will still get compared with other kernels)
    cases.push_back({{"star", 0., 0, 0}, in});
  }

  {
    Input in(SST::make_rectangle<K>(0, 0, 4, 2));
    in.with_angles(Weights{{FT(45), FT(90), FT(45), FT(90)}});
    // vertical walls are approximated with a large weight so (without snapping) the volume is not exactly 8
    // @todo add a test for snapping
    cases.push_back({{"gable roof", 8., 0, 0}, in});
  }

  return cases;
}

template <typename K>
std::vector<CGAL::Surface_mesh<typename K::Point_3> > run_kernel(const std::string& kname)
{
  using FT = typename K::FT;
  using Mesh = CGAL::Surface_mesh<typename K::Point_3>;

  using Outcome = SST::Outcome;
  using Expectations = SST::Expectations;

  std::cout << "== " << kname << " ==" << std::endl;

  std::vector<Mesh> meshes;
  for(const auto& [d, in] : kernel_cases<K>())
  {
    const std::string name = kname + ": " + d.name;
    std::cout << "  -- Test " << name << std::endl;

    meshes.emplace_back();
    Mesh& m = meshes.back();
    const Outcome o = extrude(in, m);
    CGAL_SLS_CHECK_MSG(o == Outcome::SUCCESS, name + ": " + to_string(o));

    if(o != Outcome::SUCCESS)
      continue;

    Expectations e;
    if(d.volume > 0)
      e.volume = d.volume;
    e.rel_tol = 1e-6; // vertical walls
    check_extrusion(name, in, m, e);

    if constexpr(SSE::internal::is_exact_kernel<K>)
    {
      if(d.exact_den != 0)
      {
        const FT v = CGAL::Polygon_mesh_processing::volume(m);
        CGAL_SLS_CHECK_MSG(v == FT(d.exact_num) / FT(d.exact_den), name << ": inexact volume " << v);
      }
    }
  }
  return meshes;
}

template <typename Mesh1, typename Mesh2>
void compare(const std::string& kname,
             const std::vector<Mesh1>& ms,
             const std::vector<Mesh2>& refs)
{
  const auto names = kernel_cases<EPICK>();
  for(std::size_t i=0; i<ms.size(); ++i)
  {
    if(num_faces(ms[i]) == 0 || num_faces(refs[i]) == 0)
      continue;
    const double tol = 1e-9 * SST::diag_xy(PMP::bbox(refs[i]));
    CGAL_SLS_CHECK_MSG(SST::approx_equal(PMP::volume(ms[i]), PMP::volume(refs[i]), tol),
                       kname << " vs EPECK: " << names[i].first.name);
  }
}

int main(int, char**)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  const auto ref = run_kernel<EPECK>("EPECK");
  compare("EPICK", run_kernel<EPICK>("EPICK"), ref);
  compare("EPECK_with_sqrt", run_kernel<EPECK_with_sqrt>("EPECK_with_sqrt"), ref);
  compare("SC<Exact_rational>", run_kernel<SC_rational>("Simple_cartesian<Exact_rational>"), ref);
#ifdef CGAL_USE_GMPXX
  compare("SC<mpq_class>", run_kernel<CGAL::Simple_cartesian<mpq_class> >("Simple_cartesian<mpq_class>"), ref);
#endif

  CGAL::Surface_mesh<EPECK::Point_3> m;
  CGAL_SLS_CHECK(CGAL::extrude_skeleton(SST::make_rectangle<EPICK>(0, 0, 2, 2), m));
  CGAL_SLS_CHECK(num_vertices(m) == 5);

  return SST::test_exit_code();
}
