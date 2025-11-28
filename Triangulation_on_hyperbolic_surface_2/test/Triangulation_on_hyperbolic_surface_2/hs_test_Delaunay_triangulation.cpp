#include <CGAL/Exact_rational.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Hyperbolic_surface_traits_2.h>

#include <CGAL/Delaunay_triangulation_on_hyperbolic_surface_2.h>
#include <CGAL/Triangulation_on_hyperbolic_surface_2_IO.h>

#include <iostream>
#include <sstream>
#include <vector>

using namespace CGAL;


typedef Simple_cartesian<CGAL::Exact_rational>                      Kernel;
typedef Hyperbolic_Delaunay_triangulation_traits_2<Kernel>          ParentTraits;

typedef Hyperbolic_surface_traits_2<ParentTraits>                   Traits;
typedef Hyperbolic_fundamental_domain_2<Traits>                     Domain;
typedef Delaunay_triangulation_on_hyperbolic_surface_2<Traits>      Delaunay_triangulation;
typedef typename Delaunay_triangulation::Anchor                     Anchor;
typedef typename Delaunay_triangulation::CMap                       CMap;

typedef typename Traits::FT                                         FT;
typedef typename Traits::Hyperbolic_point_2                         Point;

Domain build_domain()
{
    std::vector<Point> vertices;
    Point z0 = Point(FT("4881/5000"),FT("0"));
    Point z1 = Point(FT("9211/10000"),FT("2733/10000"));
    Point z2 = Point(FT("1709/5000"),FT("7253/10000"));
    Point z3 = Point(FT("-427262704257582473474868322141310044732400799603/1267155016747148041260345910894159385550919570000"),FT("582571804584198065321856347012850217722442509611/1267155016747148041260345910894159385550919570000"));
    Point z4 = Point(FT("-4881/5000"),FT("0"));
    Point z5 = Point(FT("-9211/10000"),FT("-2733/10000"));
    Point z6 = Point(FT("-1709/5000"),FT("-7253/10000"));
    Point z7 = Point(FT("427262704257582473474868322141310044732400799603/1267155016747148041260345910894159385550919570000"),FT("-582571804584198065321856347012850217722442509611/1267155016747148041260345910894159385550919570000"));
    vertices.push_back(z0);
    vertices.push_back(z1);
    vertices.push_back(z2);
    vertices.push_back(z3);
    vertices.push_back(z4);
    vertices.push_back(z5);
    vertices.push_back(z6);
    vertices.push_back(z7);

    std::vector<std::size_t> pairings;
    for (std::size_t k=0; k<8; ++k) {
        pairings.push_back((k+4)%8);
    }
    return Domain(vertices, pairings);
}

int main()
{
    Domain domain = build_domain();
    Delaunay_triangulation dt = Delaunay_triangulation(domain);

    assert(dt.is_valid());

    std::cout << "printing triangulation for test purposes: " << std::endl << dt;

    Anchor a0 = dt.locate(Point(0.5, 0.5)); // straight walk
    Anchor a1 = dt.locate(Point(0.5, 0.5), true); // visibility walk

    assert(dt.anchor(a0.dart).dart == dt.anchor(a1.dart).dart); //check that both anchors correspond to the same triangle

    bool same_vertices = true;
    for (unsigned i = 0; i < 3; ++i) {
        bool found = false;
        for (unsigned j = 0; j < 3; ++j) {
            if (a0.vertices[i] == a1.vertices[j]) {
                found = true;
            }
        }
        same_vertices = same_vertices && found;
    }

    assert(same_vertices);

    assert(dt.shortest_loop() != 0);
    assert(dt.shortest_non_loop_edge() != 0);

    return 0;
}
