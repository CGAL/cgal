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

    typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
    typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
    typedef Tr::Geom_traits GT;
    typedef GT::Sphere_3 Sphere_3;
    typedef GT::Point_3 Point_3;
    typedef GT::FT FT;
    typedef FT(*Function)(Point_3);
    typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;
    typedef CGAL::Surface_mesh<Point_3> Surface_mesh;
    FT sphere_function(Point_3 p) {
        const FT x2 = p.x() * p.x(), y2 = p.y() * p.y(), z2 = p.z() * p.z();
        return x2 + y2 + z2 - 1;
    }
// User-defined feature
template <typename GeomTraits, typename PointRange, typename PointMap>
class Compactness : public CGAL::Classification::Feature_base
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
    
    Image_float dtm;
    std::vector<typename FloatMap::value_type> feature_values;

public:
    Compactness(const PointRange& input, PointMap point_map, const Grid& grid, float feature_scale, const Mesh& wrap) : input(input), point_map(point_map), grid(grid)
    {
        this->set_name("Compactness");
        feature_scale *= 2.0f;

        dtm = Image_float(grid.width(), grid.height());

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


        namespace AW3 = CGAL::Alpha_wraps_3;
        typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel_Alpha;
        typedef Kernel_Alpha::Point_3 Point_Alpha;
        using Point_container = std::vector<Point_Alpha>;
        using Mesh = CGAL::Surface_mesh<Point_Alpha>;

        double start = omp_get_wtime();
        int fin = 0;
        bool print_once = false;
        //#pragma omp parallel for private(sm) shared(fin,print_once)
        #pragma omp parallel for shared(fin,print_once)
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

                    // ***************************
                    // RMD you want to replace this code for computing "Betti numbers" on the intersected mesh with code to compute the surface area and volume of intersected_mesh.
                    // ***************************



                    //compute volume 

                    double vol = 0;

                    for (Mesh::Face_range::iterator it = intersected_mesh.faces_begin(); it != intersected_mesh.faces_end(); ++it) {
                        halfedge_descriptor hd = intersected_mesh.halfedge(*it);

                        Point_container P;
                        do {
                            const Point_3& p = intersected_mesh.point(intersected_mesh.source(hd));
                            P.push_back(p);
                            hd = intersected_mesh.next(hd);
                        } while (hd != intersected_mesh.halfedge(*it));
              
            
                        // Determinant
                        double v321 = P.at(2).x() * P.at(1).y() * P.at(0).z();
                        double v231 = P.at(1).x() * P.at(2).y() * P.at(0).z();
                        double v312 = P.at(2).x() * P.at(0).y() * P.at(1).z();
                        double v132 = P.at(0).x() * P.at(2).y() * P.at(1).z();
                        double v213 = P.at(1).x() * P.at(0).y() * P.at(2).z();
                        double v123 = P.at(0).x() * P.at(1).y() * P.at(2).z();

                        double det = -v321 + v231 + v312 - v132 - v213 + v123;
                        vol += det;

                    }
                    vol /= (double)6.0f;
       
                    //compute surface area
                    double surface = 0;

                    for (Mesh::Face_range::iterator it = intersected_mesh.faces_begin(); it != intersected_mesh.faces_end(); ++it) {
                        halfedge_descriptor hd = intersected_mesh.halfedge(*it);

                        Point_container P;
                        do {
                            const Point_3& p = intersected_mesh.point(intersected_mesh.source(hd));
                            P.push_back(p);
                            hd = intersected_mesh.next(hd);
                        } while (hd != intersected_mesh.halfedge(*it));

                        double area = 0.5 * sqrt(std::pow((P.at(1).y() - P.at(0).y()) * (P.at(2).z() - P.at(0).z()) - (P.at(1).z() - P.at(0).z()) * (P.at(2).y() - P.at(0).y()), 2) +
                            (std::pow((P.at(1).z() - P.at(0).z()) * (P.at(2).x() - P.at(0).x()) - (P.at(1).x() - P.at(0).x()) * (P.at(2).z() - P.at(0).z()), 2)
                                + (std::pow((P.at(1).x() - P.at(0).x()) * (P.at(2).y() - P.at(0).y()) - (P.at(1).y() - P.at(0).y()) * (P.at(2).x() - P.at(0).x()), 2))));

                        surface += std::abs(area);
                    }
                    

                    double compactness = (3 * vol) / (feature_scale * surface);
                    if (surface == 0.f) compactness = 0.f;

                    //double compactness = 0;
                    dtm(i, j) = compactness;
                    /*
                    if (!print_once) {
                        std::cout << "(" << i << "," << j << ") " << dtm(i, j) << ", vol: " << vol << ", area: " << surface << " | " << wrap.number_of_vertices() << "," << sm.number_of_vertices() << "," << smi.number_of_vertices() << "," << intersected_mesh.number_of_vertices() << std::endl;
                        print_once = true;
                    }
                    */
                    /*
                    if (i == 0)
                        std::cout << "Surface area :" << surface << std::endl;
                    if (i == 0)
                        std::cout << "compactness :" << compactness << std::endl;

                    if (compactness == 0) {
                        CGAL::IO::write_polygon_mesh("intersected_mesh.ply", intersected_mesh, CGAL::parameters::stream_precision(17));
                    }
                    */
                }

            }

        double end = omp_get_wtime();
        printf_s("Time = %.8g\n", end - start);


        if (grid.width() * grid.height() > input.size()) {
            feature_values.resize(input.size(), 0.f);
            for (std::size_t i = 0; i < input.size(); ++i) {
                std::size_t I = grid.x(i);
                std::size_t J = grid.y(i);
                feature_values[i] = dtm(I, J);
            }
            dtm.free();
        }
       
    }
    float value(std::size_t pt_index) {
        if (feature_values.empty()) {
            std::size_t I = grid.x(pt_index);
            std::size_t J = grid.y(pt_index);
            return dtm(I, J);
        }

        return feature_values[pt_index];
    }
};
}}}