#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_3_triangulation_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>

#include <CGAL/Random.h>
#include <CGAL/point_generators_3.h>

// Traits
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Periodic_3_triangulation_traits_3<K> PK;

// Vertex type
typedef CGAL::Periodic_3_triangulation_ds_vertex_base_3<> DsVb;
typedef CGAL::Triangulation_vertex_base_3<PK,DsVb> Vb;
typedef CGAL::Alpha_shape_vertex_base_3<PK,Vb> AsVb;
// Cell type
typedef CGAL::Periodic_3_triangulation_ds_cell_base_3<> DsCb;
typedef CGAL::Triangulation_cell_base_3<PK,DsCb> Cb;
typedef CGAL::Alpha_shape_cell_base_3<PK,Cb> AsCb;

typedef CGAL::Triangulation_data_structure_3<AsVb,AsCb> Tds;
typedef CGAL::Periodic_3_Delaunay_triangulation_3<PK,Tds> P3DT3;
typedef CGAL::Alpha_shape_3<P3DT3>  Alpha_shape_3;

typedef PK::Point_3                                 Point;

int main()
{
  typedef CGAL::Creator_uniform_3<double, Point> Creator;
  CGAL::Random random(7);
  CGAL::Random_points_in_cube_3<Point, Creator> in_cube(1, random);
  std::vector<Point> pts;

  // Generating 1000 random points
  for (int i=0 ; i < 1000 ; i++) {
    Point p = *in_cube++;
    pts.push_back(p);
  }

  // Define the periodic cube
  P3DT3 pdt(PK::Iso_cuboid_3(-1,-1,-1,1,1,1));
  // Heuristic for inserting large point sets (if pts is reasonably large)
  pdt.insert(pts.begin(), pts.end(), true);
  // As pdt won't be modified anymore switch to 1-sheeted cover if possible
  if (pdt.is_triangulation_in_1_sheet()) pdt.convert_to_1_sheeted_covering();
  std::cout << "Periodic Delaunay computed." << std::endl;

  // compute alpha shape
  Alpha_shape_3 as(pdt);
  std::cout << "Alpha shape computed in REGULARIZED mode by default."
	    << std::endl;

   // find optimal alpha values
  Alpha_shape_3::NT alpha_solid = as.find_alpha_solid();
  Alpha_shape_3::Alpha_iterator opt = as.find_optimal_alpha(1);
  std::cout << "Smallest alpha value to get a solid through data points is "
	    << alpha_solid << std::endl;
  std::cout << "Optimal alpha value to get one connected component is "
	    <<  *opt    << std::endl;
  as.set_alpha(*opt);
  assert(as.number_of_solid_components() == 1);
  return 0;
}
