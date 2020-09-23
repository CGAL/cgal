#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>

typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<>          Traits;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>        PHDTriangulation;
typedef PHDTriangulation::Point                                                                                 Point;
typedef PHDTriangulation::Vertex_handle                                                                        Vertex_handle;
typedef Traits::FT                                                                                                                 FT;
typedef Traits::Circle_2                                                                                                 Circle;
typedef Traits::Euclidean_line_2                                                                                 ELine;
typedef Traits::Construct_hyperbolic_segment_2                                                         CHSegment;
typedef Traits::Construct_intersection_2                                                                 CIntersection;

int main()
{

        std::vector<Point> P(10), Pp(10);
        std::vector<ELine> S(10), Sp(10);

         P[0]  = Point(  FT(2)/FT(10),  FT(8)/FT(10));
        Pp[0] = Point( -FT(2)/FT(10), -FT(8)/FT(10));

        P[1]  = Point(  FT(4)/FT(10),  FT(8)/FT(10));
        Pp[1] = Point( -FT(4)/FT(10), -FT(8)/FT(10));

        P[2]  = Point(  FT(6)/FT(10),  FT(6)/FT(10));
        Pp[2] = Point( -FT(6)/FT(10), -FT(6)/FT(10));

        P[3]  = Point(  FT(8)/FT(10),  FT(4)/FT(10));
        Pp[3] = Point( -FT(8)/FT(10), -FT(4)/FT(10));

        P[4]  = Point(  FT(8)/FT(10),  FT(2)/FT(10));
        Pp[4] = Point( -FT(8)/FT(10), -FT(2)/FT(10));

        P[5]  = Point(  FT(8)/FT(10), -FT(2)/FT(10));
        Pp[5] = Point( -FT(8)/FT(10),  FT(2)/FT(10));

        P[6]  = Point(  FT(8)/FT(10), -FT(4)/FT(10));
        Pp[6] = Point( -FT(8)/FT(10),  FT(4)/FT(10));

        P[7]  = Point(  FT(6)/FT(10), -FT(6)/FT(10));
        Pp[7] = Point( -FT(6)/FT(10),  FT(6)/FT(10));

        P[8]  = Point(  FT(4)/FT(10), -FT(8)/FT(10));
        Pp[9] = Point( -FT(4)/FT(10),  FT(8)/FT(10));

        P[9]  = Point(  FT(2)/FT(10), -FT(8)/FT(10));
        Pp[9] = Point( -FT(2)/FT(10),  FT(8)/FT(10));

        Point O(FT(0), FT(0));

        for (int i = 0; i < 10; i++) {
                S[i]  = ELine( P[i], Pp[i] );
        }

        PHDTriangulation tri;

        Circle c1( O, FT(1)/FT(3) );
        for (int i = 0; i < 10; i++) {
                std::cout << "--- c1, i = " << i << std::endl;
                std::pair<Point, Point> res = CIntersection()(c1, S[i]);
                tri.insert(res.first);
                tri.insert(res.second);
        }

        Circle c2( O, FT(1)/FT(4) );
        for (int i = 0; i < 10; i++) {
                std::cout << "--- c2, i = " << i << std::endl;
                std::pair<Point, Point> res = CIntersection()(c2, S[i]);
                tri.insert(res.first);
                tri.insert(res.second);
        }

        Circle c3( O, FT(1)/FT(5) );
        for (int i = 0; i < 10; i++) {
                std::cout << "--- c3, i = " << i << std::endl;
                std::pair<Point, Point> res = CIntersection()(c3, S[i]);
                tri.insert(res.first);
                tri.insert(res.second);
        }

        std::cout << " -------- inserting --------" << std::endl;
        std::cout << "Vertices  in triangulation: " << tri.number_of_vertices() << std::endl;
        std::cout << "Faces     in triangulation: " << tri.number_of_faces() << std::endl;
        std::cout << "Dimension of triangulation: " << tri.dimension() << std::endl;
        assert(tri.dimension() == 2);
        assert(tri.is_valid());

        std::cout << " -------- SUCCESS --------" << std::endl;

        return 0;
}
