#ifndef CGAL_CLASSIFICATION_FEATURE_BETTI_NUMBERS_H
#define CGAL_CLASSIFICATION_FEATURE_BETTI_NUMBERS_H
#include <CGAL/Classification.h>
#include <CGAL/license/Classification.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Classification.h>
//#include <CGAL/Point_set_3.h>
//#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Real_timer.h>
#ifndef CGAL_COVERAGE_FEATURE_DISPLAY_PROGRESS
#define CGAL_COVERAGE_FEATURE_DISPLAY_PROGRESS false
#endif
// Sphere generation
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Aff_transformation_3.h>

#include <CGAL/Polygon_mesh_processing/corefinement.h>

namespace CGAL {
namespace   Classification {
namespace     Feature {

// User-defined feature
template <typename GeomTraits, typename PointRange, typename PointMap, typename Mesh, typename ConcurrencyTag>
class Betti_numbers : public CGAL::Classification::Feature_base_multi_dim
{
    using MeshKernel = typename Mesh::Point::R::Kernel;
    using MeshPoint = typename Mesh::Point;
    using Point_container = std::vector<MeshPoint>;

    // Mesh intersection
    typedef boost::graph_traits<typename Mesh>::template halfedge_descriptor        halfedge_descriptor;
    typedef boost::graph_traits<typename Mesh>::template edge_descriptor            edge_descriptor;
    typedef boost::graph_traits<typename Mesh>::template face_descriptor            face_descriptor;

    typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
    typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
    typedef Tr::Geom_traits GT;
    typedef GT::Sphere_3 Sphere_3;
    typedef GT::Point_3 Point_3;
    typedef GT::Vector_3 Vector_3;
    typedef GT::FT FT;
    typedef FT(*Function)(Point_3);
    typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;

    // Feature description
    using Grid = CGAL::Classification::Planimetric_grid<GeomTraits, PointRange, PointMap>;
    
    using Image_float = CGAL::Classification::Image<float>;

    const PointRange& input;
    PointMap point_map;
    const Grid& grid;

    Image_float dtm_0, dtm_1, dtm_2;
    std::vector<int> b0, b1, b2;


public:
    Betti_numbers(const PointRange& input, PointMap point_map, const Grid& grid, float feature_scale, const Mesh& wrap) : input(input), point_map(point_map), grid(grid)
    {
        this->set_name("Betti_numbers");
        feature_scale *= 2.0f;

        dtm_0 = Image_float(grid.width(), grid.height());
        dtm_1 = Image_float(grid.width(), grid.height());
        dtm_2 = Image_float(grid.width(), grid.height());

        // Create Sphere
        Tr tr;           // 3D-Delaunay triangulation
        C2t3 c2t3(tr);   // 2D-complex in 3D-Delaunay triangulation
         
        Surface_3 surface([] (Point_3 p) -> FT {
                              const FT x2 = p.x() * p.x(), y2 = p.y() * p.y(), z2 = p.z() * p.z();
                              return x2 + y2 + z2 - 1;
                          }, Sphere_3(CGAL::ORIGIN, 2.f)); // bounding sphere function, squared radius

        CGAL::Surface_mesh_default_criteria_3<Tr> criteria(20.f,  // angular bound  (default: 30.)
                                                           0.2,  // radius bound   (default: 0.1)
                                                           0.2); // distance bound (default: 0.1)
        // meshing surface
        CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());

        Mesh sm;
        CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, sm);

        int occupied_cells = 0;
        for (int j = 0; j < grid.height(); ++j)
            for (std::size_t i = 0; i < grid.width(); ++i) {
                if (grid.has_points(i, j)) {
                    occupied_cells++;
                }
            }

        CGAL::Cartesian_converter<GeomTraits, GT> to_tr;

        std::size_t n_cells = grid.height() * grid.width();
#if CGAL_COVERAGE_FEATURE_DISPLAY_PROGRESS
        std::cout << "No. of cells: " << grid.height() * grid.width() << " occupied cells: " << occupied_cells << std::endl;
#ifdef BOOST_TIMER_PROGRESS_DISPLAY_HPP_INCLUDED
        boost::timer::progress_display progress(n_cells);
#else
        CGAL::Real_timer t;
        t.reset();
        t.start();
#endif
#endif
        CGAL::for_each<ConcurrencyTag>
            (CGAL::make_counting_range<std::size_t>(0, n_cells),
                [&](const std::size_t& s) -> bool
                {
                    std::size_t i = s % grid.width();
                    std::size_t j = s / grid.width();
                    if (grid.has_points(i, j)) {
                        // calculate centroid of points in grid cell
                        float cx = 0, cy = 0, cz = 0;
                        PointRange::Vector_3 c;
                        int k = 0;
                        typename Grid::iterator end = grid.indices_end(i, j);
                        for (typename Grid::iterator it = grid.indices_begin(i, j); it != end; ++it) {
                            PointRange::Point_3& p = get(point_map, *(input.begin() + (*it)));
                            cx += p.x();
                            cy += p.y();
                            cz += p.z();
                            //c += PointRange::Vector_3(p.x(), p.y(), p.z());
                            k++;
                        }
                        PointRange::Point_3 centroid(cx / k, cy / k, cz / k);
                        //PointRange::Point_3 centroid(c.x() / k, c.y() / k, c.z() / k);
                        //Vector_3 centroid = to_tr(c / k);

                        Mesh smi = sm; // copy sphere mesh

                        // Translate and scale sphere mesh
                        for (Mesh::Vertex_index& p : smi.vertices()) {
                            Point_3& a = smi.point(p);
                            smi.point(p) = Point_3(centroid.x() + a.x() * feature_scale, centroid.y() + a.y() * feature_scale, centroid.z() + a.z() * feature_scale);
                        }

                        // Intersect wrap and sphere
                        Mesh intersected_mesh = wrap;
                        Mesh::Property_map<edge_descriptor, bool> is_constrained_map = intersected_mesh.add_property_map<edge_descriptor, bool>("e:is_constrained", false).first;
                        bool valid_difference = CGAL::Polygon_mesh_processing::corefine_and_compute_intersection(intersected_mesh,
                            smi,
                            intersected_mesh,
                            CGAL::parameters::default_values(),
                            CGAL::parameters::default_values(),
                            CGAL::parameters::edge_is_constrained_map(is_constrained_map));

                        // compute Betti numbers
                        int number_of_connected_components = CGAL::Polygon_mesh_processing::connected_components(intersected_mesh, intersected_mesh.add_property_map<face_descriptor, std::size_t>("f:CC").first);
                        int euler_characteristic = (int) intersected_mesh.number_of_vertices() - intersected_mesh.number_of_edges() + intersected_mesh.number_of_faces(); // assumes there are no unreferenced vertices
                        int total_genus = -euler_characteristic * .5f + number_of_connected_components; // assumes there are zero boundary cycles

                        int b_0 = number_of_connected_components;
                        int b_1 = 2 * total_genus; // normally: 2*genus + no. of boundary cycles - no. of disks
                        int b_2 = number_of_connected_components; // normally: no. of disconnected components - no. of disks (assumes no voids)

                        dtm_0(i, j) = (float)b_0;
                        dtm_1(i, j) = (float)b_1;
                        dtm_2(i, j) = (float)b_2;
                    }
#if CGAL_COVERAGE_FEATURE_DISPLAY_PROGRESS
#ifdef BOOST_TIMER_PROGRESS_DISPLAY_HPP_INCLUDED
                    ++progress;
#endif
#endif
                    return true;
                });

#if CGAL_COVERAGE_FEATURE_DISPLAY_PROGRESS
#ifndef BOOST_TIMER_PROGRESS_DISPLAY_HPP_INCLUDED
                t.stop();
                std::cout << "Took " << t.time() << " s." << std::endl;
#endif
#endif

        if (grid.width() * grid.height() > input.size()) {
            b0.resize(input.size(), 0.f);
            b1.resize(input.size(), 0.f);
            b2.resize(input.size(), 0.f);
            for (std::size_t i = 0; i < input.size(); ++i) {
                std::size_t I = grid.x(i);
                std::size_t J = grid.y(i);
                b0[i] = (int) dtm_0(I, J);
                b1[i] = (int) dtm_1(I, J);
                b2[i] = (int) dtm_2(I, J);
            }
        }
    }

    float value(std::size_t pt_index, std::size_t dim_index)
    {
        if (b0.empty() && b1.empty() && b2.empty()) {
            std::size_t I = grid.x(pt_index);
            std::size_t J = grid.y(pt_index);
            if (dim_index == 0)
                return dtm_0(I, J);
            else if (dim_index == 1)
                return dtm_1(I, J);
            else if (dim_index == 2)
                return dtm_2(I, J);
        } else {
            if (dim_index == 0)
                return b0.at(pt_index);
            else if (dim_index == 1)
                return b1.at(pt_index);
            else if (dim_index == 2)
                return b2.at(pt_index);
        }
        return 0;
    }
};

} // namespace Feature
} // namespace Classification
} // namespace CGAL
#endif CGAL_CLASSIFICATION_FEATURE_BETTI_NUMBERS_H