#include <CGAL/Hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/point_generators_2.h>

typedef CGAL::Hyperbolic_Delaunay_triangulation_traits_2<>  Traits;
typedef CGAL::Hyperbolic_Delaunay_triangulation_2<Traits>        HDTriangulation;
typedef HDTriangulation::Point                                                                 Point;
typedef HDTriangulation::Vertex_handle                                                 Vertex_handle;

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
        std::vector<Vertex_handle> verts;
        for (unsigned int i = 0; i < pts.size(); ++i)
                verts.push_back(tri.insert(pts[i]));

        std::size_t nv = tri.number_of_vertices();
        std::size_t nf = tri.number_of_hyperbolic_faces();

        std::cout << " -------- inserting --------" << std::endl;
        std::cout << "Vertices  in triangulation: " << nv << std::endl;
        std::cout << "Faces     in triangulation: " << nf << std::endl;
        std::cout << "Dimension of triangulation: " << tri.dimension() << std::endl;
        assert(tri.dimension() == 2);

        std::cout << " -------- removing --------" << std::endl;
        for (unsigned int i = 0; i < verts.size(); ++i)
        {
                tri.remove(verts[i]);
        }

        std::cout << "Vertices  in triangulation: " << tri.number_of_vertices() << std::endl;
        std::cout << "Faces     in triangulation: " << tri.number_of_hyperbolic_faces()    << std::endl;
        std::cout << "Dimension of triangulation: " << tri.dimension() << std::endl;

        assert(tri.is_valid());

        std::cout << " -------- SUCCESS --------" << std::endl;

        HDTriangulation tri2;
        verts.clear();
        for (unsigned int i = 0; i < pts.size(); ++i)
                verts.push_back(tri2.insert(pts[i]));
        std::cout << " -------- inserting --------" << std::endl;
        std::cout << "Vertices  in triangulation: " << tri2.number_of_vertices()                  << std::endl;
        std::cout << "Faces     in triangulation: " << tri2.number_of_hyperbolic_faces() << std::endl;
        std::cout << "Dimension of triangulation: " << tri2.dimension()                                  << std::endl;
        assert(tri2.dimension() == 2);

        tri2.remove(verts.begin(), verts.end());

        std::cout << "Vertices  in triangulation: " << tri2.number_of_vertices() << std::endl;
        std::cout << "Faces     in triangulation: " << tri2.number_of_hyperbolic_faces()    << std::endl;
        std::cout << "Dimension of triangulation: " << tri2.dimension() << std::endl;

        assert(tri2.is_valid());

        std::cout << " -------- SUCCESS --------" << std::endl;


        return 0;
}
