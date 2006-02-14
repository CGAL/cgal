// file: examples/Fixed_precision_nt/delaunay.C

#include <CGAL/Cartesian.h>

#include <cstdio>
#include <iostream>

#include <CGAL/Triangulation_short_names_2.h>

#include <CGAL/Fixed_precision_nt.h>
#include <CGAL/Timer.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/point_generators_2.h>

typedef CGAL::Fixed_precision_nt coord_type;

static bool Fixed_precision_nt_init_result =
            CGAL::Fixed_precision_nt::init(2000.0);

typedef CGAL::Cartesian<coord_type>               Repclass;
typedef Repclass::Point_2                         Point_;
typedef CGAL::Delaunay_triangulation_2<Repclass>  Delaunay_;

int main(int argc, char* argv[])
{
    CGAL::force_ieee_double_precision();

    CGAL::Timer t;
    Delaunay_ D;

    CGAL::Random_points_in_square_2<Point_,
                                    CGAL::Creator_uniform_2<double,Point_> >
        Input ( 1.0 );

    int N;
    if (argc==2)
      CGAL_CLIB_STD::sscanf(argv[1], "%d", &N); 
    else {
      N=100;
      std::cerr<<"usage : "<<argv[0]<<" nb-of-points"<<std::endl<<std::endl;
    }
    std::cout << "Delaunay of "<<N<<" random points"<<std::endl;
    t.start();
    while(N--)
      D.insert( *++Input);
    t.stop();
    std::cout << " in "<<t.time()<<" seconds"<<std::endl;
    return 0;
}
