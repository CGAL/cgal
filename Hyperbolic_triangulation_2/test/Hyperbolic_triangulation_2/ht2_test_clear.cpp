#include <CGAL/Hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/point_generators_2.h>

typedef CGAL::Hyperbolic_Delaunay_triangulation_traits_2<>  Traits;
typedef CGAL::Hyperbolic_Delaunay_triangulation_2<Traits>        HDTriangulation;
typedef HDTriangulation::Point                                                                 Point;

typedef CGAL::Random_points_in_square_2<
            Point,
             CGAL::Creator_uniform_2< double, Point > >                         Point_generator;

int main()
{

        Point_generator gen(1.);
        std::vector<Point> pts;
        for (int i = 0; i < 100; ++i)
                pts.push_back(*(++gen));

        HDTriangulation tri;
        tri.insert(pts.begin(), pts.end());

        std::size_t nv = tri.number_of_vertices();
        std::size_t nf = tri.number_of_hyperbolic_faces();

        std::cout << " -------- inserting --------" << std::endl;
        std::cout << "Vertices  in triangulation: " << nv << std::endl;
        std::cout << "Faces     in triangulation: " << nf << std::endl;
        std::cout << "Dimension of triangulation: " << tri.dimension() << std::endl;
        assert(tri.dimension() == 2);

        tri.clear();
        std::cout << " -------- clearing --------" << std::endl;
        std::cout << "Vertices  in triangulation: " << tri.number_of_vertices() << std::endl;
        std::cout << "Faces     in triangulation: " << tri.number_of_hyperbolic_faces()    << std::endl;
        std::cout << "Dimension of triangulation: " << tri.dimension() << std::endl;
        assert(tri.number_of_vertices() == 0);
        assert(tri.number_of_hyperbolic_faces() == 0);
        assert(tri.dimension() < 2);

        tri.insert(pts.begin(), pts.end());
        std::cout << " -------- re-inserting --------" << std::endl;
        std::cout << "Vertices  in triangulation: " << tri.number_of_vertices() << std::endl;
        std::cout << "Faces     in triangulation: " << tri.number_of_hyperbolic_faces()    << std::endl;
        std::cout << "Dimension of triangulation: " << tri.dimension() << std::endl;
        assert(tri.number_of_hyperbolic_faces() == nf);
        assert(tri.number_of_vertices() == nv);
        assert(tri.dimension() == 2);

        assert(tri.is_valid());

        std::cout << " -------- SUCCESS --------" << std::endl;

        return 0;
}
