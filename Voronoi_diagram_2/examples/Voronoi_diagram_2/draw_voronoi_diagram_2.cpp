// standard includes
#include <fstream>

// includes for drawing the Voronoi Diagram
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>
#include <CGAL/draw_voronoi_diagram_2.h>

// typedefs for defining the adaptor
typedef CGAL::Exact_predicates_inexact_constructions_kernel                  K;
typedef CGAL::Delaunay_triangulation_2<K>                                    DT;
typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT>                 AT;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT> AP;
typedef CGAL::Voronoi_diagram_2<DT,AT,AP>                                    VD;

// typedef for the result type of the point location
typedef AT::Site_2                    Site_2;

int main(int argc, char* argv[])
{
  VD vd;
  std::ifstream ifs((argc>1)?argv[1]:"data/data4.dt.cin");
  assert(ifs);

  Site_2 t;
  while ( ifs >> t ) { vd.insert(t); }
  ifs.close();

  assert( vd.is_valid() );

  CGAL::draw(vd);

  return EXIT_SUCCESS;
}
