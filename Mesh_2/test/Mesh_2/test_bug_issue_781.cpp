#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

// mesh refinement
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_2                                     Point;

typedef CGAL::Triangulation_vertex_base_2<Kernel>           Vb;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel> Fbb;
typedef CGAL::Delaunay_mesh_face_base_2<Kernel,Fbb>         Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>         TDS;
typedef CGAL::Exact_intersections_tag                       Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel,
						   TDS,
						   Itag>    CDT;

typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>            Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria>              Mesher;

int main(int, char* [])
{
  std::cerr.precision(17);
  CDT cdt;

  Point pt(5.4691594172333904, 44.256641611715409);

  cdt.insert(pt);

  Point s(5.4688567999999904, 44.256729199999995);
  Point t(5.4693788499999929, 44.256578099999999);

  cdt.insert_constraint(s, t);

  assert(cdt.is_valid());
  Mesher mesher(cdt);
  mesher.set_criteria(Criteria(0.125));

  std::cout << "refine mesh..." << std::flush;
  mesher.refine_mesh();
  std::cout << " complete: "
	    << cdt.number_of_vertices() << " vertices" << std::endl;
  assert(cdt.number_of_vertices() == 3);
  return 0;
}
