#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Periodic_3_regular_triangulation_traits_3.h>
#include <CGAL/Periodic_3_regular_triangulation_3.h>

#include <CGAL/Random.h>

#include <fstream>
#include <iostream>

// Traits
typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Periodic_3_regular_triangulation_traits_3<K>   PK;

// Vertex type
typedef CGAL::Periodic_3_triangulation_ds_vertex_base_3<>    DsVb;
typedef CGAL::Regular_triangulation_vertex_base_3<PK,DsVb>   Vb;
typedef CGAL::Alpha_shape_vertex_base_3<PK,Vb>               AsVb;
// Cell type
typedef CGAL::Periodic_3_triangulation_ds_cell_base_3<>      DsCb;
typedef CGAL::Regular_triangulation_cell_base_3<PK,DsCb>     Cb;
typedef CGAL::Alpha_shape_cell_base_3<PK,Cb>                 AsCb;

typedef CGAL::Triangulation_data_structure_3<AsVb,AsCb>      Tds;
typedef CGAL::Periodic_3_regular_triangulation_3<PK,Tds>     P3RT3;
typedef CGAL::Alpha_shape_3<P3RT3>                           Alpha_shape_3;

typedef P3RT3::Bare_point                                    Bare_point;
typedef P3RT3::Weighted_point                                Weighted_point;

int main()
{
  CGAL::Random random(8);
  std::list<Weighted_point> pts;

  // read input
  std::ifstream is("./data/bunny_1000");
  int n;
  is >> n;
  std::cout << "Reading " << n << " points " << std::endl;
  for( ; n>0 ; n--)
  {
    Bare_point bp;
    if(is >> bp)
      pts.emplace_back(bp, 0.0001 * random.get_double(0., 0.015625)); // arbitrary weights
  }

  // Define the periodic cube
  P3RT3 prt(PK::Iso_cuboid_3(-0.1,0.,-0.1, 0.1,0.2,0.1));

  // Heuristic for inserting large point sets (if pts is reasonably large)
  prt.insert(pts.begin(), pts.end(), true);

  // As prt won't be modified anymore switch to 1-sheeted cover if possible
  if(prt.is_triangulation_in_1_sheet())
  {
    std::cout << "Switching to 1-sheeted covering" << std::endl;
    prt.convert_to_1_sheeted_covering();
  }

  std::cout << "Periodic Regular computed." << std::endl;

  // compute alpha shape
  Alpha_shape_3 as(prt);
  std::cout << "Alpha shape computed in REGULARIZED mode by default." << std::endl;

   // find optimal alpha values
  Alpha_shape_3::NT alpha_solid = as.find_alpha_solid();
  Alpha_shape_3::Alpha_iterator opt = as.find_optimal_alpha(1);
  std::cout << "Smallest alpha value to get a solid through data points is " << alpha_solid << std::endl;
  std::cout << "Optimal alpha value to get one connected component is " << *opt << std::endl;
  as.set_alpha(*opt);
  assert(as.number_of_solid_components() == 1);

  return 0;
}
