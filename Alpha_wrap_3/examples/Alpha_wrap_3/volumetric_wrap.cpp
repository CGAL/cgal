#include "output_helper.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/alpha_wrap_3.h>

#include <CGAL/tetrahedral_remeshing.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_cell_base_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_vertex_base_3.h>
#include <CGAL/Simplicial_mesh_cell_base_3.h>
#include <CGAL/Simplicial_mesh_vertex_base_3.h>

#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/property_map.h>
#include <CGAL/Real_timer.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/draw_triangulation_3.h>

#include <CGAL/IO/Triangulation_off_ostream_3.h>
#include <CGAL/IO/File_medit.h>

#include <iostream>
#include <string>

namespace PMP = CGAL::Polygon_mesh_processing;
namespace AW3i = CGAL::Alpha_wraps_3::internal;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;

using Points = std::vector<Point_3>;
using Face = std::array<std::size_t, 3>;
using Faces = std::vector<Face>;

using Mesh = CGAL::Surface_mesh<Point_3>;

// If we provide a triangulation, AW3 uses its Gt, so we have to make the Gt stack explicit
using Gtb = AW3i::Alpha_wrap_AABB_geom_traits<K>; // provides Ball_3
using Gt = CGAL::Robust_circumcenter_filtered_traits_3<Gtb>; // better inexact constructions (not mandatory)

// Since we are going to use tetrahedral remeshing on the underlying triangulation,
// we need special vertex and cell base types that meets the requirements of the
// tetrahedral remeshing concepts
using Vbbb = AW3i::Alpha_wrap_triangulation_vertex_base_3<K>;
using Vbb = CGAL::Simplicial_mesh_vertex_base_3<K, int, int, int, int, Vbbb>;
using Vb = CGAL::Tetrahedral_remeshing::Remeshing_vertex_base_3<K, Vbb>;

using Cbbb = AW3i::Alpha_wrap_triangulation_cell_base_3<K>;
using Cbb = CGAL::Simplicial_mesh_cell_base_3<K, int, int, Cbbb>;
using Cb = CGAL::Tetrahedral_remeshing::Remeshing_cell_base_3<K, Cbb>;

using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;

using Delaunay_triangulation = CGAL::Delaunay_triangulation_3<Gt, Tds, CGAL::Fast_location>;

// because the Fast_location does all kinds of rebinding shenanigans + T3_hierarchy is in the stack...
using Triangulation = CGAL::Triangulation_3<typename Delaunay_triangulation::Geom_traits,
                                            typename Delaunay_triangulation::Triangulation_data_structure>;

using Facet = Triangulation::Facet;

int main(int argc, char** argv)
{
  // Read the input
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/bull.off");
  std::cout << "Reading " << filename << "..." << std::endl;

  Points points;
  Faces faces;
  if(!CGAL::IO::read_polygon_soup(filename, points, faces) || faces.empty())
  {
    std::cerr << "Invalid input:" << filename << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Input: " << points.size() << " vertices, " << faces.size() << " faces" << std::endl;

  // Compute the alpha and offset values
  const double relative_alpha = (argc > 2) ? std::stod(argv[2]) : 20.;
  const double relative_offset = (argc > 3) ? std::stod(argv[3]) : 600.;

  CGAL::Bbox_3 bbox;
  for(const Point_3& p : points)
    bbox += p.bbox();

  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()) +
                                       CGAL::square(bbox.zmax() - bbox.zmin()));

  const double alpha = diag_length / relative_alpha;
  const double offset = diag_length / relative_offset;
  std::cout << "alpha: " << alpha << ", offset: " << offset << std::endl;

  // Construct the wrap
  CGAL::Real_timer t;
  t.start();

  using Oracle = CGAL::Alpha_wraps_3::internal::Triangle_soup_oracle<K>;

  Oracle oracle(K{});
  oracle.add_triangle_soup(points, faces, CGAL::parameters::default_values());

  CGAL::Alpha_wraps_3::internal::Alpha_wrapper_3<Oracle, Delaunay_triangulation> aw3(oracle);
  Mesh wrap;
  aw3(alpha, offset, wrap);

  t.stop();
  std::cout << "Result: " << num_vertices(wrap) << " vertices, " << num_faces(wrap) << " faces" << std::endl;
  std::cout << "Took " << t.time() << " s." << std::endl;

  // Get the interior tetrahedrization
  auto dt = aw3.triangulation();

  // Save the result
  const std::string output_name = generate_output_name(filename, relative_alpha, relative_offset);
  std::cout << "Writing to " << output_name << std::endl;
  CGAL::IO::write_polygon_mesh(output_name, wrap, CGAL::parameters::stream_precision(17));

  // Remesh the interior of the wrap
  const Delaunay_triangulation& aw3_dt = aw3.triangulation();
  const Triangulation& aw3_tr = static_cast<const Triangulation&>(aw3_dt);
  Triangulation tr = aw3_tr; // intentional copy

  std::cout << "BEFORE: " << tr.number_of_vertices() << " vertices, " << tr.number_of_cells() << " cells" << std::endl;

  // Set up the c3t3 information
  for(auto v : tr.finite_vertex_handles())
    v->set_dimension(3);

  for(auto c : tr.finite_cell_handles())
  {
    if(c->is_outside())
      c->set_subdomain_index(0);
    else
      c->set_subdomain_index(1);

    // if the neighboring cell has a different outside info, put the vertices
    // of the common face on the surface boundary
    for(int i=0; i<4; ++i)
    {
      if(c->neighbor(i)->is_outside() != c->is_outside())
      {
        c->set_surface_patch_index(i, 1);
        for(int j=1; j<4; ++j)
          c->vertex((i+j)%4)->set_dimension(2);
      }
    }
  }

  std::ofstream out_before("before_remeshing.mesh");
  CGAL::IO::write_MEDIT(out_before, tr);

  // edge length of equilateral triangle with circumradius alpha
  // const double l = 2 * alpha * 0.8660254037844386; // sqrt(3)/2

  // edge length of regular tetrahedron with circumradius alpha
  const double l = 1.6329931618554521 * alpha; // sqrt(8/3)

  CGAL::tetrahedral_isotropic_remeshing(tr, l,
                                        CGAL::parameters::remesh_boundaries(false)
                                                         .number_of_iterations(5));

  std::cout << "AFTER: " << tr.number_of_vertices() << " vertices, " << tr.number_of_cells() << " cells" << std::endl;

  std::ofstream out_after("after_remeshing.mesh");
  CGAL::IO::write_MEDIT(out_after, tr);

  return EXIT_SUCCESS;
}
