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

#include <CGAL/IO/Triangulation_off_ostream_3.h>
#include <CGAL/IO/File_medit.h>


#include <Mesh_smoothing_3/Mesh_smoothing_3.h>

#include <iostream>
#include <string>
#include <chrono>

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

using SC = CGAL::Simple_cartesian<double>;
using SC_Point_3 = SC::Point_3;
using SC_Vector_3 = SC::Vector_3;
using SC_Iso_cuboid_3 = SC::Iso_cuboid_3;
using SC2GT = CGAL::Cartesian_converter<SC, Delaunay_triangulation::Geom_traits>;

class Tetrahedral_mesh_wrapper {
public:
    using Cell_descriptor = Triangulation::Cell_handle;
    using Vertex_descriptor =  Triangulation::Vertex_handle;
    using Point_3 = K::Point_3;

    std::size_t nb_cells() const { return tetmesh.number_of_finite_cells(); }
    std::size_t nb_vertices() const { return tetmesh.number_of_vertices(); }

    Point_3 vertex_coordinates(Vertex_descriptor vertex) const { return tetmesh.point(vertex); }
    void set_new_vertex_coordinates(Vertex_descriptor vertex, Point_3 coord) {
        vertex->set_point(coord);
    }  // only non const

    auto cell_range() const {
        return cell_vector_range;
    }
    auto cell_vertices(Cell_descriptor cell) const {
        std::array<Vertex_descriptor, 4> vertices;
        for (int i = 0; i < 4; ++i) {
            vertices[static_cast<unsigned>(i)] = cell->vertex(i);
        }
        return vertices;
    }
    std::array<Point_3, 4> cell_reference_shape(Cell_descriptor) const { return Mesh_smoothing_3::Shapes::VTK_TETRAHEDRON<Point_3>(); }
public:
    Tetrahedral_mesh_wrapper(Triangulation &tetmesh_, std::set<int> regions)
    : tetmesh(tetmesh_)
    {
        cell_vector_range.reserve(tetmesh.number_of_cells());
        for(auto c : tetmesh.finite_cell_handles()) {
            if (regions.contains(c->subdomain_index())) {
                cell_vector_range.push_back(c);
            }
        }

    }
    Triangulation &tetmesh;
    std::vector<Triangulation::Cell_handle> cell_vector_range;
};

class Triangle_boundary_wrapper {
public:
    using Face_descriptor = unsigned;
    using Vertex_descriptor = Triangulation::Vertex_handle;
    using Normal_3 =  K::Vector_3;
    using Surface_patch_index = unsigned;
    std::size_t nb_faces() const { return faces.size(); }
    auto face_range() const {
        return Mesh_smoothing_3::utils::Contiguous_unsigned_range{0, faces.size()};
    }
    unsigned patch_id(Face_descriptor) const { return 0; }
    std::size_t nb_face_vertices(Face_descriptor) const { return 3; }
    auto face_vertices(Face_descriptor face) const {
        return faces[face];
    }

public:
    std::vector<std::array<Vertex_descriptor, 3>> const &faces;
};



int main(int argc, char** argv)
{
    // Read the input
    const std::string filename = (argc > 1) ? argv[1] : "../../data/joint.off";
    std::cout << "Reading " << filename << "..." << std::endl;

    Points points;
    Faces faces;
    if(!CGAL::IO::read_polygon_soup(filename, points, faces, CGAL::parameters::verbose(true)) || faces.empty())
    {
        std::cerr << "Invalid input:" << filename << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Input: " << points.size() << " vertices, " << faces.size() << " faces" << std::endl;

    // Compute the alpha and offset values
    const double relative_alpha = (argc > 2) ? std::stod(argv[2]) : 60.;
    const double relative_offset = (argc > 3) ? std::stod(argv[3]) : 600.;

    CGAL::Bbox_3 bbox;
    for(const Point_3& p : points) bbox += p.bbox();

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



    Delaunay_triangulation& aw3_dt = aw3.triangulation();
    for(auto c : aw3_dt.finite_cell_handles())
    {
        if(c->is_outside())
            c->set_subdomain_index(0);
        else
            c->set_subdomain_index(1);
    }

    const Triangulation& aw3_tr = static_cast<const Triangulation&>(aw3_dt);
    Triangulation tr = aw3_tr; // intentional copy

    std::cout << "BEFORE REMESHING: " << tr.number_of_vertices() << " vertices, " << tr.number_of_cells() << " cells" << std::endl;

    // Set up the c3t3 information
    for(auto v : tr.finite_vertex_handles()) v->set_dimension(3);

    for(auto c : tr.finite_cell_handles())
    {
        for(int i=0; i<4; ++i)
        {
            if(c->neighbor(i)->subdomain_index() != c->subdomain_index())
            {
                c->set_surface_patch_index(i, 1);
                for(int j=1; j<4; ++j) {
                    c->vertex((i+j)%4)->set_dimension(2);
                }
            }
        }
    }

    // edge length of regular tetrahedron with circumradius alpha
    const double l = 1.6329931618554521 * alpha; // sqrt(8/3)

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    // remeshing the interior
    CGAL::tetrahedral_isotropic_remeshing(tr, l, CGAL::parameters::remesh_boundaries(false).number_of_iterations(5));  
    // Remeshing the surface may lead to worse quality around features because of point equidistribution on the surface. 
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count() << "[ms]" << std::endl;
    std::cout << "AFTER REMESHING: " << tr.number_of_vertices() << " vertices, " << tr.number_of_cells() << " cells" << std::endl;

    auto boundary_query = [&](Point_3 pt, unsigned, double) -> std::tuple<Point_3, K::Vector_3, double>{
        auto pp_and_prim = oracle.tree().closest_point_and_primitive(pt);
        Point_3 proj = pp_and_prim.first;
        K::Vector_3 normal { proj, pt };
        normal = normal / CGAL::approximate_sqrt(normal * normal);
        proj = proj + offset * normal;
        return {proj, normal, 1.};
    };


    std::vector<std::array<Triangulation::Vertex_handle, 3>> boundary_faces;
    std::unordered_map<Triangulation::Vertex_handle, std::array<bool, 3>> locked;
    std::vector<Triangulation::Vertex_handle> boundary_vertices;
    for(auto v : tr.finite_vertex_handles()) {
        locked.emplace(v, std::array<bool, 3>{true, true, true});
    }
    for(auto c : tr.finite_cell_handles())
    {
        if (c->subdomain_index() != 1) continue;

        for(int i=0; i<4; ++i)
        {
            if(c->neighbor(i)->subdomain_index() != 1)
            {
                for(int j=0; j<4; ++j) {
                    locked.at(c->vertex(j)) = {false, false, false};
                }
                boundary_faces.push_back({});
                for(int j=1; j<4; ++j) {
                    boundary_faces.back()[j-1] = c->vertex((i+j)%4);
                    boundary_vertices.push_back(c->vertex((i+j)%4));
                }
            }

        }
    } 

    Tetrahedral_mesh_wrapper mesh_wrapper(tr, {1,2});
    Triangle_boundary_wrapper boundary_wrapper {boundary_faces};

    Mesh_smoothing_3::Mesh_smoother smoother(mesh_wrapper, boundary_wrapper);
    smoother.set_verbose();

    smoother.set_locked_vertices(boundary_vertices);
    smoother.run();

    std::ofstream all_after_smoothed("offset_wrapper.mesh");
    CGAL::IO::write_MEDIT(all_after_smoothed, tr, CGAL::parameters::all_cells(false));

    smoother.clear_locks();
    smoother.set_vertices_dim_locks(locked);
    smoother.set_boundary_query(boundary_query);
    smoother.set_max_number_of_iteration(100);
    smoother.run();
    smoother.clear_locks();
    smoother.run();


    std::ofstream inside_after_smoothed("offset_ours.mesh");
    CGAL::IO::write_MEDIT(inside_after_smoothed, tr, CGAL::parameters::all_cells(false));
    return EXIT_SUCCESS;
}
