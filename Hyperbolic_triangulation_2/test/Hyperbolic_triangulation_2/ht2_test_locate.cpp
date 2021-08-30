#include <CGAL/Hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/point_generators_2.h>

typedef CGAL::Hyperbolic_Delaunay_triangulation_traits_2<>  Traits;
typedef CGAL::Hyperbolic_Delaunay_triangulation_2<Traits>        HDTriangulation;
typedef HDTriangulation::FT                                                                 FT;
typedef HDTriangulation::Point                                                                 Point;
typedef HDTriangulation::Face_handle                                                 Face_handle;
typedef HDTriangulation::Locate_type                                                 Locate_type;
typedef HDTriangulation::Hyperbolic_segment                                 Hyperbolic_segment;

int main()
{
        FT F100(100);
        Point p1(-FT(81)/F100, -FT(35)/F100        );
        Point p2(-FT(78)/F100,  FT(12)/F100        );
        Point p3(-FT(64)/F100,  FT(59)/F100        );
        Point p4(-FT(31)/F100, -FT(25)/F100        );
        Point p5( FT(66)/F100,  FT(18)/F100 );

        HDTriangulation tri;
        tri.insert(p1);
        tri.insert(p2);
        tri.insert(p3);
        tri.insert(p4);
        tri.insert(p5);

        std::size_t nv = tri.number_of_vertices();
        std::size_t nf = tri.number_of_hyperbolic_faces();

        std::cout << " -------- inserting --------" << std::endl;
        std::cout << "Vertices  in triangulation: " << nv << std::endl;
        std::cout << "Faces     in triangulation: " << nf << std::endl;
        std::cout << "Dimension of triangulation: " << tri.dimension() << std::endl;
        assert(tri.dimension() == 2);

    int idx;
    Face_handle fh;

        // Test a vertex location
        Locate_type lt1;
    fh = tri.locate(p1, lt1, idx);
        assert(lt1 == HDTriangulation::VERTEX);

        // Test an edge location
        Hyperbolic_segment s1 = tri.geom_traits().construct_hyperbolic_segment_2_object()(p1, p2);
        Hyperbolic_segment b1 = tri.geom_traits().construct_hyperbolic_bisector_2_object()(p1, p2);
        Point i1 = tri.geom_traits().construct_intersection_2_object()(s1, b1);
        Locate_type lt2;
        fh = tri.locate(i1, lt2, idx);
        assert(lt2 == HDTriangulation::EDGE);

        // Test a "dangling" edge location
        Hyperbolic_segment s2 = tri.geom_traits().construct_hyperbolic_segment_2_object()(p4, p5);
        Hyperbolic_segment b2 = tri.geom_traits().construct_hyperbolic_bisector_2_object()(p4, p5);
        Point i2 = tri.geom_traits().construct_intersection_2_object()(s2, b2);
        Locate_type lt3;
        fh = tri.locate(i2, lt3, idx);
        assert(lt3 == HDTriangulation::EDGE);

        std::cout << " -------- SUCCESS --------" << std::endl;

        return 0;
}
