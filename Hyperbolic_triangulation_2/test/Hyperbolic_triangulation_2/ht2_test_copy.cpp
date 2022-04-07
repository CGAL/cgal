#include <CGAL/Hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/point_generators_2.h>

typedef CGAL::Hyperbolic_Delaunay_triangulation_traits_2<>  Traits;
typedef CGAL::Hyperbolic_Delaunay_triangulation_2<Traits>        HDTriangulation;
typedef HDTriangulation::Point                                                                 Point;
typedef HDTriangulation::Vertex_handle                                                 Vertex_handle;

typedef CGAL::Random_points_in_disc_2<
            Point,
             CGAL::Creator_uniform_2< double, Point > >                         Point_generator;

int main()
{

        Point_generator gen(0.95);
        std::vector<Point> pts1;
        for (int i = 0; i < 75; ++i)
                pts1.push_back(*(++gen));

        HDTriangulation tri1;
        tri1.insert(pts1.begin(), pts1.end());
        std::cout << "Vertices  in tri1: " << tri1.number_of_vertices()                 << std::endl;
        std::cout << "Faces     in tri1: " << tri1.number_of_hyperbolic_faces() << std::endl;

        HDTriangulation tri2;
        tri2 = tri1;
        std::cout << "Vertices  in tri2: " << tri2.number_of_vertices()                 << std::endl;
        std::cout << "Faces     in tri2: " << tri2.number_of_hyperbolic_faces() << std::endl;

        assert(tri1.is_valid());
        assert(tri2.is_valid());
        assert(tri1.number_of_vertices() == tri2.number_of_vertices());
        assert(tri1.number_of_hyperbolic_faces() == tri2.number_of_hyperbolic_faces());

        std::cout << " -------- SUCCESS --------" << std::endl;


        return 0;
}
