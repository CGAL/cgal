#include <CGAL/Classification.h>
#include <CGAL/license/Classification.h>

#include <omp.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Classification.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Real_timer.h>
// Sphere generation
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Aff_transformation_3.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/alpha_wrap_3.h>

#include <CGAL/Polygon_mesh_processing/corefinement.h>

namespace CGAL {
namespace   Classification {
namespace     Feature {

class Betti_0 : public CGAL::Classification::Feature_base {
    std::vector<int> b0;
public:
    Betti_0(std::vector<int>& b0) : b0(b0) {
        this->set_name("betti0");
    }

    float value(std::size_t pt_index) {
        return (float) b0[pt_index];
    }
};

class Betti_1 : public CGAL::Classification::Feature_base {
    std::vector<int> b1;
public:
    Betti_1(std::vector<int>& b1) : b1(b1) {
        this->set_name("betti1");
    }

    float value(std::size_t pt_index) {
        return (float) b1[pt_index];
    }
};

class Betti_2 : public CGAL::Classification::Feature_base {
    std::vector<int> b2;
public:
    Betti_2(std::vector<int>& b2) : b2(b2) {
        this->set_name("betti2");
    }

    float value(std::size_t pt_index) {
        return (float) b2[pt_index];
    }

};

typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;
typedef FT(*Function)(Point_3);
typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;
FT sphere_function_betti(Point_3 p) {
    const FT x2 = p.x() * p.x(), y2 = p.y() * p.y(), z2 = p.z() * p.z();
    return x2 + y2 + z2 - 1;
}

template <typename GeomTraits, typename PointRange, typename PointMap, typename Grid>
std::tuple< std::vector<int>, std::vector<int>, std::vector<int> > compute_betti_numbers(const PointRange& input, PointMap point_map, const Grid& grid, float feature_scale, const Mesh& wrap){//(const PointRange& input, PointMap point_map, const float feature_scale, const Mesh& wrap) {
    // Mesh intersection
    typedef boost::graph_traits<Mesh>::halfedge_descriptor        halfedge_descriptor;
    typedef boost::graph_traits<Mesh>::edge_descriptor            edge_descriptor;
    typedef boost::graph_traits<Mesh>::face_descriptor            face_descriptor;

    using Image_float = CGAL::Classification::Image<float>;

    using FloatMap = typename PointRange::template Property_map<float>;

    Image_float dtm_0, dtm_1, dtm_2;
    //std::vector<typename FloatMap::value_type> b0, b1, b2;
    std::vector<int> b0, b1, b2;



    feature_scale *= 2.0f;

    dtm_0 = Image_float(grid.width(), grid.height());
    dtm_1 = Image_float(grid.width(), grid.height());
    dtm_2 = Image_float(grid.width(), grid.height());

    // Create Sphere
    Tr tr;            // 3D-Delaunay triangulation
    C2t3 c2t3(tr);   // 2D-complex in 3D-Delaunay triangulation

    // defining the surface
    Surface_3 surface(sphere_function_betti,             // pointer to function
        Sphere_3(CGAL::ORIGIN, 2.)); // bounding sphere
    // Note that "2." above is the *squared* radius of the bounding sphere!
    // defining meshing criteria

    CGAL::Surface_mesh_default_criteria_3<Tr> criteria(20.,  // angular bound  (default: 30.)
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
    std::cout << "No. of cells: " << grid.height() * grid.width() << " occupied cells: " << occupied_cells << std::endl;
    //boost::timer::progress_display show_progress(grid.height() * grid.width());

    double start = omp_get_wtime();
    int fin = 0;
#pragma omp parallel for private(sm) shared(fin)
    for (int k = 0; k < grid.height() * grid.width(); ++k) {
        std::size_t i = k % grid.width();
        std::size_t j = k / grid.width();
        //for (int j = 0; j < grid.height(); ++j)
        //    for (std::size_t i = 0; i < grid.width(); ++i) {
        if (grid.has_points(i, j)) {
            if (++fin % int(occupied_cells / 100.f) == 0) {
                //std::cout << float(fin) / occupied_cells * 100.f << "%" << std::endl;
                printf("\r%.2f%%", float(fin) / occupied_cells * 100.f);
            }
            // calculate centroid of points in grid cell
            float cx = 0, cy = 0, cz = 0;
            //Vector_3 c;
            int k = 0;
            typename Grid::iterator end = grid.indices_end(i, j);
            for (typename Grid::iterator it = grid.indices_begin(i, j); it != end; ++it) {
                PointRange::Point_3& p = get(point_map, *(input.begin() + (*it)));

                cx += p.x();
                cy += p.y();
                cz += p.z();
                //c += p;
                k++;
            }
            PointRange::Point_3 centroid(cx / k, cy / k, cz / k);
            //PointRange::Point_3 centroid(c.x() / k, c.y() / k, c.z() / k);

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
                CGAL::parameters::edge_is_constrained_map(intersected_mesh.property_map<edge_descriptor, bool>("e:is_constrained").first));

            // compute Betti numbers
            int number_of_connected_components = CGAL::Polygon_mesh_processing::connected_components(intersected_mesh, intersected_mesh.add_property_map<face_descriptor, std::size_t>("f:CC").first);
            int euler_characteristic = (int)intersected_mesh.number_of_vertices() - intersected_mesh.number_of_edges() + intersected_mesh.number_of_faces(); // assumes there are no unreferenced vertices
            int total_genus = -euler_characteristic * .5f + number_of_connected_components; // assumes there are zero boundary cycles

            int b_0 = number_of_connected_components;
            int b_1 = 2 * total_genus; // normally: 2*genus + no. of boundary cycles - no. of disks
            int b_2 = number_of_connected_components; // normally: no. of disconnected components - no. of disks

            dtm_0(i, j) = (float)b_0;
            dtm_1(i, j) = (float)b_1;
            dtm_2(i, j) = (float)b_2;
        }

    }

    double end = omp_get_wtime();
    printf_s("Time = %.8g\n", end - start);

    //if (grid.width() * grid.height() > input.size()) {
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

    dtm_0.free();
    dtm_1.free();
    dtm_2.free();
    //}

    return std::make_tuple(b0, b1, b2);
}
/*
// User-defined feature
template <typename GeomTraits, typename PointRange, typename PointMap>
class BettiNumbers : public CGAL::Classification::Feature_base
{
    

    // Mesh intersection
    typedef boost::graph_traits<Mesh>::halfedge_descriptor        halfedge_descriptor;
    typedef boost::graph_traits<Mesh>::edge_descriptor            edge_descriptor;
    typedef boost::graph_traits<Mesh>::face_descriptor            face_descriptor;
    //namespace PMP = CGAL::Polygon_mesh_processing;
    //namespace params = CGAL::parameters;



    using Grid = CGAL::Classification::Planimetric_grid<GeomTraits, PointRange, PointMap>;
    
    using Image_float = CGAL::Classification::Image<float>;

    using FloatMap = typename PointRange::template Property_map<float>;

    const PointRange& input;
    PointMap point_map;
    const Grid& grid;
    
    Image_float dtm_0, dtm_1, dtm_2;
    std::vector<typename FloatMap::value_type> b0, b1, b2;

public:
    std::tuple< std::vector<int>, std::vector<int>, std::vector<int> > ComputeBettiNumbers(const PointRange& input, PointMap point_map, const Grid& grid, float feature_scale, const Mesh& wrap) : input(input), point_map(point_map), grid(grid)
    {
        this->set_name("Betti_numbers");
        feature_scale *= 2.0f;

        dtm_0 = Image_float(grid.width(), grid.height());
        dtm_1 = Image_float(grid.width(), grid.height());
        dtm_2 = Image_float(grid.width(), grid.height());

        // Create Sphere
        Tr tr;            // 3D-Delaunay triangulation
        C2t3 c2t3(tr);   // 2D-complex in 3D-Delaunay triangulation

        // defining the surface
        Surface_3 surface(sphere_function,             // pointer to function
                          Sphere_3(CGAL::ORIGIN, 2.)); // bounding sphere
        // Note that "2." above is the *squared* radius of the bounding sphere!
        // defining meshing criteria

        CGAL::Surface_mesh_default_criteria_3<Tr> criteria(20.,  // angular bound  (default: 30.)
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
        std::cout << "No. of cells: " << grid.height() * grid.width() << " occupied cells: " << occupied_cells << std::endl;
        //boost::timer::progress_display show_progress(grid.height() * grid.width());

        double start = omp_get_wtime();
        int fin = 0;
#pragma omp parallel for private(sm) shared(fin)
        for (int k = 0; k < grid.height() * grid.width(); ++k) {
            std::size_t i = k % grid.width();
            std::size_t j = k / grid.width();
            //for (int j = 0; j < grid.height(); ++j)
            //    for (std::size_t i = 0; i < grid.width(); ++i) {
            if (grid.has_points(i, j)) {
                if (++fin % int(occupied_cells / 100.f) == 0) {
                    //std::cout << float(fin) / occupied_cells * 100.f << "%" << std::endl;
                    printf("%.2f%%\r", float(fin) / occupied_cells * 100.f);
                }
                // calculate centroid of points in grid cell
                float cx = 0, cy = 0, cz = 0;
                //Vector_3 c;
                int k = 0;
                typename Grid::iterator end = grid.indices_end(i, j);
                for (typename Grid::iterator it = grid.indices_begin(i, j); it != end; ++it) {
                    PointRange::Point_3& p = get(point_map, *(input.begin() + (*it)));

                    cx += p.x();
                    cy += p.y();
                    cz += p.z();
                    //c += p;
                    k++;
                }
                PointRange::Point_3 centroid(cx / k, cy / k, cz / k);
                //PointRange::Point_3 centroid(c.x() / k, c.y() / k, c.z() / k);

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
                    CGAL::parameters::edge_is_constrained_map(intersected_mesh.property_map<edge_descriptor, bool>("e:is_constrained").first));

                // compute Betti numbers
                int number_of_connected_components = CGAL::Polygon_mesh_processing::connected_components(mi, mi.add_property_map<face_descriptor, std::size_t>("f:CC").first);
                int euler_characteristic = (int)mi.number_of_vertices() - mi.number_of_edges() + mi.number_of_faces(); // assumes there are no unreferenced vertices
                int total_genus = -euler_characteristic * .5f + number_of_connected_components; // assumes there are zero boundary cycles

                int b_0 = number_of_connected_components;
                int b_1 = 2 * total_genus; // normally: 2*genus + no. of boundary cycles - no. of disks
                int b_2 = number_of_connected_components; // normally: no. of disconnected components - no. of disks

                dtm_0(i, j) = (float)b_0;
                dtm_1(i, j) = (float)b_1;
                dtm_2(i, j) = (float)b_2;
            }

        }

        double end = omp_get_wtime();
        printf_s("Time = %.8g\n", end - start);

        //if (grid.width() * grid.height() > input.size()) {
        b0.resize(input.size(), 0.f);
        b1.resize(input.size(), 0.f);
        b2.resize(input.size(), 0.f);
        for (std::size_t i = 0; i < input.size(); ++i) {
            std::size_t I = grid.x(i);
            std::size_t J = grid.y(i);
            feature_values[i] = dtm(I, J);
            b0[i] = dtm_0(I, J);
            b1[i] = dtm_1(I, J);
            b2[i] = dtm_2(I, J);
        }

        dtm.free();
        //}

        return std::make_tuple(b0, b1, b2);
    }
    float value(std::size_t pt_index) {
        return -1; // RMD you should not be using this method
    }
};*/


}}}