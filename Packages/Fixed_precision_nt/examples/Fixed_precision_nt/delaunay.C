#include <CGAL/basic.h>
#include <stdio.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>
#include <strstream.h>

#include <CGAL/Triangulation_short_names_2.h>

#include <CGAL/Fixed_precision_nt.h>
#include <CGAL/Timer.h>

#include <CGAL/Cartesian.h>
#include <CGAL/squared_distance_2.h>   // to avoid a g++ problem
#include <CGAL/Point_2.h>

#include <CGAL/Triangulation_euclidean_traits_2.h>

#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/point_generators_2.h>

using namespace CGAL;

typedef Fixed_precision_nt coord_type;
static bool Fixed_precision_nt_init_result
             = Fixed_precision_nt::init(2000.0);

typedef Cartesian<coord_type>  Repclass;

typedef Point_2<Repclass>  Point_;

typedef Triangulation_euclidean_traits_2<Repclass> Traits_;
typedef Triangulation_vertex_base_2<Traits_> Vb;
typedef Triangulation_face_base_2<Traits_>  Fb;
typedef Triangulation_default_data_structure_2<Traits_,Vb,Fb> Tds;
typedef Delaunay_triangulation_2<Traits_,Tds>  Delaunay_;


int main(int argc, char* argv[])
{
    Timer t;
    Delaunay_ D;

    Random_points_in_square_2<Point_,
      Creator_uniform_2<double,Point_> > Input ( 1.0 );

    int N; if (argc==2) sscanf(argv[1], "%d", &N); 
    else {N=100; cerr<<"usage : "<<argv[0]<<" nb-of-points"<<endl<<endl;}
    cout << "Delaunay of "<<N<<" random points"<<endl;
    t.start();
    while(N--) D.insert( *++Input);
    t.stop();
    cout << " in "<<t.time()<<" seconds"<<endl;
    return 0;
}

