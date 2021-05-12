#include <CGAL/Hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/point_generators_2.h>

typedef CGAL::Hyperbolic_Delaunay_triangulation_traits_2<>  Traits;
typedef CGAL::Hyperbolic_Delaunay_triangulation_2<Traits>        HDTriangulation;
typedef HDTriangulation::FT                                                                 FT;
typedef HDTriangulation::Point                                                                 Point;
typedef HDTriangulation::Vertex_handle                                                 Vertex_handle;
typedef HDTriangulation::Vertex_circulator                                         Vertex_circulator;


int main()
{
        HDTriangulation tri;

        std::cout << "Constructing points..." << std::endl;
        FT F100(100);
        Point p1 = tri.geom_traits().construct_point_2_object()(-FT(62)/F100, -FT(49)/F100);
        Point p2 = tri.geom_traits().construct_point_2_object()(-FT(57)/F100,  FT(68)/F100);
        Point p3 = tri.geom_traits().construct_point_2_object()( FT(60)/F100,  FT(58)/F100);
        Point p4 = tri.geom_traits().construct_point_2_object()( FT(48)/F100, -FT(64)/F100);
        Point p5 = tri.geom_traits().construct_point_2_object()(-FT(07)/F100, -FT(04)/F100);

        std::cout << "Inserting points..." << std::endl;
        Vertex_handle v1 = tri.insert(p1);
        Vertex_handle v2 = tri.insert(p2);
        Vertex_handle v3 = tri.insert(p3);
        Vertex_handle v4 = tri.insert(p4);
        Vertex_handle v5 = tri.insert(p5);

        std::cout << "Testing first circulator..." << std::endl;
        Vertex_circulator vc1 = tri.adjacent_vertices(v5);

        assert(vc1 == v1 || vc1 == v2 || vc1 == v3 || vc1 == v4);
        ++vc1;
        assert(vc1 == v1 || vc1 == v2 || vc1 == v3 || vc1 == v4);
        ++vc1;
        assert(vc1 == v1 || vc1 == v2 || vc1 == v3 || vc1 == v4);
        ++vc1;
        assert(vc1 == v1 || vc1 == v2 || vc1 == v3 || vc1 == v4);
        --vc1;
        assert(vc1 == v1 || vc1 == v2 || vc1 == v3 || vc1 == v4);
        --vc1;
        assert(vc1 == v1 || vc1 == v2 || vc1 == v3 || vc1 == v4);
        --vc1;
        assert(vc1 == v1 || vc1 == v2 || vc1 == v3 || vc1 == v4);
        --vc1;
        assert(vc1 == v1 || vc1 == v2 || vc1 == v3 || vc1 == v4);
        std::cout << " -------- SUCCESS --------" << std::endl;


        std::cout << std::endl << "Testing second circulator..." << std::endl;
        Vertex_circulator vc2 = tri.adjacent_vertices(v1);
        assert(vc2 == v5);
        ++vc2;
        assert(vc2 == v5);
        --vc2;
        assert(vc2 == v5);

        std::cout << " -------- SUCCESS --------" << std::endl;

        return 0;
}
