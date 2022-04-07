#include <iostream>
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>
#include <CGAL/algorithm.h>
#include <CGAL/double.h>
#include <CGAL/Simple_cartesian.h>
#include <cassert>
typedef double RT;

typedef CGAL::Simple_cartesian<RT>                                K;
typedef K::Point_3                                                 Point_3;
typedef K::Vector_3                                                 Vector_3;
typedef K::Tetrahedron_3                                         Tetrahedron_3;
typedef std::vector<Point_3>                                         Container;
typedef CGAL::Random_points_in_tetrahedron_3<Point_3>                 Point_generator;

template<class InputIterator>
bool inside_or_close_to_tetrahedron(const Tetrahedron_3& tet,InputIterator begin, InputIterator end) {
        while(begin!=end) {
                Tetrahedron_3 OABC = Tetrahedron_3(*begin, tet.vertex(0),
                                tet.vertex(1), tet.vertex(2));
                Tetrahedron_3 OABD = Tetrahedron_3(*begin, tet.vertex(0),
                                tet.vertex(1), tet.vertex(3));
                Tetrahedron_3 OBCD = Tetrahedron_3(*begin, tet.vertex(1),
                                tet.vertex(2), tet.vertex(3));
                Tetrahedron_3 OACD = Tetrahedron_3(*begin, tet.vertex(0),
                                tet.vertex(2), tet.vertex(3));
                K::FT OABC_volume = fabs(OABC.volume());
                K::FT OABD_volume = fabs(OABD.volume());
                K::FT OBCD_volume = fabs(OBCD.volume());
                K::FT OACD_volume = fabs(OACD.volume());
                K::FT tet_volume = fabs(tet.volume());
                if
                        (fabs(OABC_volume+OABD_volume+OBCD_volume+OACD_volume-tet_volume)>1e-15)
                {
                        return false;
                }
                ++begin;
        }
        return true;
}

int main() {
        CGAL::Random rand;
        Container point_set;
        const int number_tetrahedrons = 10;
        const int number_points = 1000;
        for(int i = 0; i < number_tetrahedrons; ++i) {
                Point_3 pts[4];
                for(int j = 0; j < 4; ++j) {
                        pts[j]=Point_3(rand.get_double(),rand.get_double(),rand.get_double());
                }
                Tetrahedron_3 tet(pts[0],pts[1],pts[2],pts[3]);
                Point_generator g1( pts[0], pts[1], pts[2], pts[3] );
                Point_generator g2( tet );
                Point_generator g3( g1 );

                point_set.clear();
                std::copy_n( g1, number_points,
                               std::back_inserter(point_set));
                assert(inside_or_close_to_tetrahedron(tet,point_set.begin(),point_set.end()));

                point_set.clear();
                std::copy_n( g2, number_points,
                               std::back_inserter(point_set));
                assert(inside_or_close_to_tetrahedron(tet,point_set.begin(),point_set.end()));

                point_set.clear();
                std::copy_n( g3, number_points,
                               std::back_inserter(point_set));
                assert(inside_or_close_to_tetrahedron(tet,point_set.begin(),point_set.end()));
        }
        return 0;
}
