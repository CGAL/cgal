#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/interpolation_functions.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Generalized_barycentric_coordinates_2.h>

// Some convenient typedefs.

// General.
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Projection_traits_xy_3<Kernel>                Projection;

typedef Projection::FT      Scalar;
typedef Projection::Point_2 Point;

typedef std::vector<Scalar> Scalar_vector;
typedef std::vector<Point>  Point_vector;

// Coordinates related.
typedef CGAL::Barycentric_coordinates::Mean_value_2<Projection>                                      Mean_value;
typedef CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Mean_value, Projection> Mean_value_coordinates;

// Triangulation related.
typedef CGAL::Delaunay_mesh_face_base_2<Projection>                  Face_base;
typedef CGAL::Triangulation_vertex_base_2<Projection>                Vertex_base;
typedef CGAL::Triangulation_data_structure_2<Vertex_base, Face_base> TDS;
typedef CGAL::Constrained_Delaunay_triangulation_2<Projection, TDS>  CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>                     Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria>                       Mesher;
typedef CDT::Finite_vertices_iterator                                Vertex_iterator;
typedef CDT::Vertex_handle                                           Vertex_handle;

// Interpolation related.
typedef CGAL::Interpolation_traits_2<Projection>                             Interpolation_traits;
typedef CGAL::Data_access< std::map<Point, Scalar, Projection::Less_xy_2 > > Value_access;

// STD.
using std::cout; using std::endl; using std::string;

int main()
{
    // Construct a polygon bounding a piece of three-dimensional terrain.
    // Note that z-coordinate of each vertex represents the height function.
    // Projection in 2D is done automatically by the Projection traits class.
    const int number_of_vertices = 50;
    Point_vector vertices(number_of_vertices);

    vertices[0]  = Point(0.03, 0.05, 0.000); vertices[1]  = Point(0.07, 0.04, 10.00); vertices[2]  = Point(0.10, 0.04, 20.00);
    vertices[3]  = Point(0.14, 0.04, 30.00); vertices[4]  = Point(0.17, 0.07, 40.00); vertices[5]  = Point(0.19, 0.09, 50.00);
    vertices[6]  = Point(0.22, 0.11, 60.00); vertices[7]  = Point(0.25, 0.11, 70.00); vertices[8]  = Point(0.27, 0.10, 80.00);
    vertices[9]  = Point(0.30, 0.07, 90.00); vertices[10] = Point(0.31, 0.04, 100.0); vertices[11] = Point(0.34, 0.03, 110.0);
    vertices[12] = Point(0.37, 0.02, 120.0); vertices[13] = Point(0.40, 0.03, 130.0); vertices[14] = Point(0.42, 0.04, 140.0);
    vertices[15] = Point(0.44, 0.07, 150.0); vertices[16] = Point(0.45, 0.10, 160.0); vertices[17] = Point(0.46, 0.13, 170.0);
    vertices[18] = Point(0.46, 0.19, 180.0); vertices[19] = Point(0.47, 0.26, 190.0); vertices[20] = Point(0.47, 0.31, 200.0);
    vertices[21] = Point(0.47, 0.35, 210.0); vertices[22] = Point(0.45, 0.37, 220.0); vertices[23] = Point(0.41, 0.38, 230.0);
    vertices[24] = Point(0.38, 0.37, 240.0); vertices[25] = Point(0.35, 0.36, 250.0); vertices[26] = Point(0.32, 0.35, 260.0);
    vertices[27] = Point(0.30, 0.37, 270.0); vertices[28] = Point(0.28, 0.39, 280.0); vertices[29] = Point(0.25, 0.40, 290.0);
    vertices[30] = Point(0.23, 0.39, 300.0); vertices[31] = Point(0.21, 0.37, 310.0); vertices[32] = Point(0.21, 0.34, 320.0);
    vertices[33] = Point(0.23, 0.32, 330.0); vertices[34] = Point(0.24, 0.29, 340.0); vertices[35] = Point(0.27, 0.24, 350.0);
    vertices[36] = Point(0.29, 0.21, 360.0); vertices[37] = Point(0.29, 0.18, 370.0); vertices[38] = Point(0.26, 0.16, 380.0);
    vertices[39] = Point(0.24, 0.17, 390.0); vertices[40] = Point(0.23, 0.19, 400.0); vertices[41] = Point(0.24, 0.22, 410.0);
    vertices[42] = Point(0.24, 0.25, 420.0); vertices[43] = Point(0.21, 0.26, 430.0); vertices[44] = Point(0.17, 0.26, 440.0);
    vertices[45] = Point(0.12, 0.24, 450.0); vertices[46] = Point(0.07, 0.20, 460.0); vertices[47] = Point(0.03, 0.15, 470.0);
    vertices[48] = Point(0.01, 0.10, 480.0); vertices[49] = Point(0.02, 0.07, 490.0);

    // Mesh this polygon.

    // Create a constrained Delaunay triangulation.  
    CDT cdt;

    std::vector<Vertex_handle> vertex_handle(number_of_vertices);

    // Insert vertices of the polygon as our initial point set.
    for(int i = 0; i < number_of_vertices; ++i) vertex_handle[i] = cdt.insert(vertices[i]);

    // Insert constraints - edges of the polygon - in order to mesh only the polygon's interior.
    for(int i = 0; i < number_of_vertices; ++i) cdt.insert_constraint(vertex_handle[i], vertex_handle[(i + 1) % number_of_vertices]);

    Mesher mesher(cdt);

    // Set a criteria on how to mesh.
    mesher.set_criteria(Criteria(0.01, 0.01));

    // Mesh the polygon.
    mesher.refine_mesh();

    // Compute mean value coordinates and use them to interpolate data from the polygon's boundary to its interior.

    // Associate each point with the corresponding function value and coordinates.
    std::map<Point, Scalar, Projection::Less_xy_2> point_function_value;
    std::vector< std::pair<Point, Scalar> >        point_coordinates(number_of_vertices);

    for(int i = 0; i < number_of_vertices; ++i)
            point_function_value.insert(std::make_pair(vertices[i], vertices[i].z()));

    // Create an instance of the class with mean value coordinates.
    Mean_value_coordinates mean_value_coordinates(vertices.begin(), vertices.end());

    // Store all new interior points with interpolated data here.
    std::vector<Point> points(cdt.number_of_vertices());

    cout << endl << "Result of the height interpolation: " << endl << endl;

    // Compute coordinates and interpolate the boundary data to the polygon's interior.
    int index = 0;
    for(Vertex_iterator vertex_iterator = cdt.finite_vertices_begin(); vertex_iterator != cdt.finite_vertices_end(); ++vertex_iterator) {
        Scalar_vector coordinates;

        const Point &point = vertex_iterator->point();
        mean_value_coordinates(point, std::back_inserter(coordinates));

        for(int j = 0; j < number_of_vertices; ++j) 
            point_coordinates[j] = std::make_pair(vertices[j], coordinates[j]);

        Scalar f = CGAL::linear_interpolation(point_coordinates.begin(), point_coordinates.end(), Scalar(1), Value_access(point_function_value));
        points[index] = Point(point.x(), point.y(), f);
        cout << "The interpolated height with index " << index << " is " << f << ";" << endl;
        ++index;
    }
    cout << endl;

    return EXIT_SUCCESS;
}
