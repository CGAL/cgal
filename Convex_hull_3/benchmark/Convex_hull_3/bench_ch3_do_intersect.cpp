#define CGAL_PROFILE_CONVEX_HULL_DO_INTERSECT

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Extreme_points_traits_adapter_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/convex_hull_with_hierarchy.h>
#include <CGAL/Convex_hull_3/predicates.h>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <CGAL/optimal_bounding_box.h>

#include <CGAL/Random.h>
#include <CGAL/Timer.h>
#include <CGAL/Real_timer.h>

#include <vector>
#include <fstream>

template <typename K>
struct Test
{
  typedef typename K::FT              FT;
  typedef typename K::Point_3         P;
  typedef typename K::Segment_3       S;
  typedef typename K::Vector_3        V;
  typedef typename K::Ray_3           R;
  typedef typename K::Line_3          L;
  typedef typename K::Triangle_3      T;
  typedef typename K::Plane_3         Pl;
  typedef typename K::Tetrahedron_3   Tet;
  typedef typename K::Iso_cuboid_3    Cub;
  typedef typename std::vector<P>     PR;

private:
  CGAL::Random& r;
  double m = 0, M = 1;
  const double PI=3.14159265358979323846;

public:
   Test(CGAL::Random& r) : r(r) { }

private:
  P random_point() const
  {
    return P(FT(r.get_double(m, M)), FT(r.get_double(m, M)), FT(r.get_double(m, M)));
  }

  P random_sphere_point() const
  {
    FT a=r.get_double(0,2*PI);
    FT b=r.get_double(-PI/2,PI/2);
    return P( std::cos(a)*std::cos(b), std::sin(a)*std::cos(b), std::sin(b));
  }

  void Tet_tet(int N)
  {
    std::cout << "Tetrahedron likely to intersect" << std::endl;
    std::cout << "Do intersect of package Intersection_3" << std::endl;
    CGAL::Real_timer t;
    t.start();
    for(int i=0; i<N; ++i)
    {
      P p0 = random_point();
      P p1 = random_point();
      P p2 = random_point();
      P p3 = random_point();

      P q0 = random_point();
      P q1 = random_point();
      P q2 = random_point();
      P q3 = random_point();

      CGAL::do_intersect(Tet(p0, p1, p2, p3), Tet(q0, q1, q2, q3));
    }
    t.stop();
    std::cout << t.time() << " sec" << std::endl;

    std::cout << "Do intersect of Convex Hull" << std::endl;
    t.reset();
    t.start();
    for(int i=0; i<N; ++i)
    {
      PR a({random_point(), random_point(), random_point(), random_point()});
      PR b({random_point(), random_point(), random_point(), random_point()});

      CGAL::Convex_hull_3::do_intersect(a, b);
    }
    t.stop();
    std::cout << t.time() << " sec" << std::endl << std::endl;

  }

  void Tet_mirror(int N)
  {
    std::cout << "Tetrahedrons closed" << std::endl;
    std::cout << "Do intersect of package Intersection_3" << std::endl;
    CGAL::Real_timer t;
    t.start();
    for(int i=0; i<N; ++i)
    {
      P p0 = random_point();
      if(p0.x() > p0.y())
        p0 = P(-p0.x(), -p0.y(), -p0.z());
      P p1 = random_point();
      if(p1.x() > p1.y())
        p1 = P(-p1.x(), -p1.y(), -p1.z());
      P p2 = random_point();
      if(p2.x() > p2.y())
        p2 = P(-p2.x(), -p2.y(), -p2.z());
      P p3 = random_point();
      if(p3.x() > p3.y())
        p3 = P(-p3.x(), -p3.y(), -p3.z());

      P q0 = random_point();
      if(q0.x() < q0.y())
        q0 = P(-q0.x(), -q0.y(), -q0.z());
      P q1 = random_point();
      if(q1.x() < q1.y())
        q1 = P(-q1.x(), -q1.y(), -q1.z());
      P q2 = random_point();
      if(q2.x() < q2.y())
        q2 = P(-q2.x(), -q2.y(), -q2.z());
      P q3 = random_point();
      if(q3.x() < q3.y())
        q3 = P(-q3.x(), -q3.y(), -q3.z());

      CGAL::do_intersect(Tet(p0, p1, p2, p3), Tet(q0, q1, q2, q3));
    }
    t.stop();
    std::cout << t.time() << " sec" << std::endl;

    std::cout << "Do intersect of Convex Hull" << std::endl;
    t.reset();
    t.start();
    for(int i=0; i<N; ++i)
    {
      P p0 = random_point();
      if(p0.x() > p0.y())
        p0 = P(-p0.x(), -p0.y(), -p0.z());
      P p1 = random_point();
      if(p1.x() > p1.y())
        p1 = P(-p1.x(), -p1.y(), -p1.z());
      P p2 = random_point();
      if(p2.x() > p2.y())
        p2 = P(-p2.x(), -p2.y(), -p2.z());
      P p3 = random_point();
      if(p3.x() > p3.y())
        p3 = P(-p3.x(), -p3.y(), -p3.z());

      P q0 = random_point();
      if(q0.x() < q0.y())
        q0 = P(-q0.x(), -q0.y(), -q0.z());
      P q1 = random_point();
      if(q1.x() < q1.y())
        q1 = P(-q1.x(), -q1.y(), -q1.z());
      P q2 = random_point();
      if(q2.x() < q2.y())
        q2 = P(-q2.x(), -q2.y(), -q2.z());
      P q3 = random_point();
      if(q3.x() < q3.y())
        q3 = P(-q3.x(), -q3.y(), -q3.z());

      PR a({p0,p1,p2,p3});
      PR b({q0,q1,q2,q3});

      CGAL::Convex_hull_3::do_intersect(a, b);
    }
    t.stop();
    std::cout << t.time() << " sec" << std::endl << std::endl;

  }

  void Tet_stretched(int N)
  {
    std::cout << "Tetrahedrons streched" << std::endl;
    std::cout << "Do intersect of package Intersection_3" << std::endl;
    CGAL::Real_timer t;
    t.start();
    for(int i=0; i<N; ++i)
    {
      P p0 = random_point()+V(-1,1,0);
      P p1 = random_point()+V(-1,1,0);
      P p2 = random_point()+V(1,-1,0);
      P p3 = random_point()+V(1,-1,0);

      P q0 = random_point()+V(-1,-1,0);
      P q1 = random_point()+V(-1,-1,0);
      P q2 = random_point()+V(1,1,0);
      P q3 = random_point()+V(1,1,0);

      CGAL::do_intersect(Tet(p0, p1, p2, p3), Tet(q0, q1, q2, q3));
    }
    t.stop();
    std::cout << t.time() << " sec" << std::endl;

    std::cout << "Do intersect of Convex Hull" << std::endl;
    t.reset();
    t.start();
    for(int i=0; i<N; ++i)
    {
      PR a({random_point()+V(-1,1,0),  random_point()+V(-1,1,0),  random_point()+V(1,-1,0), random_point()+V(1,-1,0)});
      PR b({random_point()+V(-1,-1,0), random_point()+V(-1,-1,0), random_point()+V(1,1,0),  random_point()+V(1,1,0)});

      CGAL::Convex_hull_3::do_intersect(a, b);
    }
    t.stop();
    std::cout << t.time() << " sec" << std::endl << std::endl;

  }

  void Tet_shift(int N)
  {
    std::cout << "Tetrahedrons that do not intersect" << std::endl;
    std::cout << "Do intersect of package Intersection_3" << std::endl;
    CGAL::Real_timer t;
    t.start();
    for(int i=0; i<N; ++i)
    {
      P p0 = random_point();
      P p1 = random_point();
      P p2 = random_point();
      P p3 = random_point();

      P q0 = random_point()+V(1,1,1);
      P q1 = random_point()+V(1,1,1);
      P q2 = random_point()+V(1,1,1);
      P q3 = random_point()+V(1,1,1);

      CGAL::do_intersect(Tet(p0, p1, p2, p3), Tet(q0, q1, q2, q3));
    }
    t.stop();
    std::cout << t.time() << " sec" << std::endl;

    std::cout << "Do intersect of Convex Hull" << std::endl;
    t.reset();
    t.start();
    for(int i=0; i<N; ++i)
    {
      PR a({random_point(), random_point(), random_point(), random_point()});
      PR b({random_point()+V(1,1,1), random_point()+V(1,1,1), random_point()+V(1,1,1), random_point()+V(1,1,1)});

      CGAL::Convex_hull_3::do_intersect(a, b);
    }
    t.stop();
    std::cout << t.time() << " sec\n" << std::endl;

  }

  void Test_sphere(int N, int M)
  {
    std::cout << "Sphere of size " << M << std::endl;
    CGAL::Real_timer t;
    std::vector< std::vector<P> > points_ranges;
    std::vector< CGAL::Surface_mesh<P> > hulls;
    std::vector< CGAL::Convex_hull_with_hierarchy<P> > hulls_wth_hierarchy;

    double nb_tests=double(N*(N-1)/2);

    for(int i=0; i<N; ++i)
    {
      std::vector<P> sp;
      V vec( P(0,0,0), random_sphere_point());
      vec*=1.2;
      for(int i=0; i<M; ++i)
        sp.push_back(random_sphere_point()+vec);
      points_ranges.push_back(std::move(sp));
    }

    std::cout << "Compute convex hulls" << std::endl;
    t.start();
    for(int i=0; i<N; ++i)
    {
      CGAL::Surface_mesh<P> hull;
      CGAL::convex_hull_3(points_ranges[i].begin(), points_ranges[i].end(), hull);
      hulls.push_back(std::move(hull));
    }
    t.stop();
    std::cout << "  " << t.time() << " sec" << std::endl;
    t.reset();

    std::cout << "Comstruct hierarchy of convex hulls" << std::endl;
    t.start();
    for(int i=0; i<N; ++i)
      hulls_wth_hierarchy.emplace_back(hulls[i]);
    t.stop();
    std::cout << "  " << t.time() << " sec" << std::endl;
    t.reset();

    std::cout << "Do intersect with RangePoints" << std::endl;
    nb_visited=0;
    t.start();
    for(int i=0; i<N; ++i)
      for(int j=i+1; j<N; ++j)
        CGAL::Convex_hull_3::do_intersect(points_ranges[i], points_ranges[j]);
    t.stop();
    std::cout << "  Number points read per test: " << double(nb_visited)/nb_tests << std::endl;
    std::cout << "  " << t.time()/nb_tests << " sec per test" << std::endl;
    t.reset();

    std::cout << "Do intersect with Hulls" << std::endl;
    nb_visited=0;
    t.start();
    for(int i=0; i<N; ++i)
      for(int j=i+1; j<N; ++j)
        CGAL::Convex_hull_3::do_intersect(hulls[i], hulls[j]);
    t.stop();
    std::cout << "  Number points read per test: " << double(nb_visited)/nb_tests << std::endl;
    std::cout << "  " << t.time()/nb_tests << " sec per test" << std::endl;
    t.reset();

    std::cout << "Do intersect with Hierarchical Hulls" << std::endl;
    nb_visited=0;
    t.start();
    for(int i=0; i<N; ++i)
      for(int j=i+1; j<N; ++j)
        CGAL::Convex_hull_3::do_intersect(hulls_wth_hierarchy[i], hulls_wth_hierarchy[j]);
    t.stop();

    std::cout << "  Number points read per test: ";
    std::cout << double(nb_visited)/nb_tests << " ";
    std::cout << std::endl;
    std::cout << "  " << t.time()/nb_tests << " sec per test" << std::endl;
    t.reset();
    std::cout << std::endl;
  }

public:
  void run()
  {
    std::cout << "Kernel: " << typeid(K).name() << std::endl;
    Tet_tet(1000);
    Tet_mirror(1000);
    Tet_shift(1000);
    Tet_stretched(1000);

    Test_sphere(1000, 40);
    Test_sphere(1000, 100);
    Test_sphere(1000, 350);
    Test_sphere(1000, 1000);
    Test_sphere(300, 3500);
    Test_sphere(50, 10000);
    Test_sphere(25, 35000);
    Test_sphere(25, 100000);
    Test_sphere(10, 350000);
    Test_sphere(10, 1000000);
    std::cout << std::endl;
  }
};

template <typename K>
void bench_on_data(const std::string &f1, const std::string &f2){
  typedef typename K::Point_3                                               Point_3;
  typedef CGAL::Surface_mesh<Point_3>                              Mesh;

  std::vector<typename K::Point_3> pts1, pts2;
  std::vector<boost::container::small_vector<std::size_t, 3>> trs1, trs2;
  if (!CGAL::IO::read_polygon_soup(f1, pts1, trs1))
  {
    std::cerr << "Cannot read " << f1 << "\n";
  }
  if (!CGAL::IO::read_polygon_soup(f2, pts2, trs2))
  {
    std::cerr << "Cannot read " << f2 << "\n";
  }

  CGAL::Real_timer t;
  t.start();
  Mesh hull1, hull2;
  CGAL::convex_hull_3(pts1.begin(), pts1.end(), hull1);
  CGAL::convex_hull_3(pts2.begin(), pts2.end(), hull2);
  t.stop();
  std::cout << "Computing convex hulls: " << t.time() << " sec" << std::endl;
  std::cout << "Convex hull size:" << vertices(hull1).size() << ", " << vertices(hull2).size() << "\n" << std::endl;
  t.reset();

  t.start();
  CGAL::Convex_hull_3::do_intersect(hull1, hull2);
  t.stop();
  std::cout << "Do intersect: " << t.time() << " sec\n" << std::endl;
  t.reset();

  t.start();
  std::array<Point_3, 8> obb1, obb2;
  CGAL::oriented_bounding_box(hull1, obb1, CGAL::parameters::use_convex_hull(false));
  CGAL::oriented_bounding_box(hull2, obb2, CGAL::parameters::use_convex_hull(false));
  t.stop();
  std::cout << "Computing Obbs: " << t.time() << " sec\n" << std::endl;
  t.reset();

  t.start();
  CGAL::Convex_hull_3::do_intersect(obb1, obb2);
  t.stop();
  std::cout << "Do intersect with Obbs: " << t.time() << " sec\n" << std::endl;
}

int main(int argc, char** argv)
{
  // std::cout.precision(17);
  // std::cerr.precision(17);

  std::cout << "3D Distance tests" << std::endl;

  CGAL::Random rp;
  CGAL::Random r(argc==1?rp.get_seed():std::stoi(argv[1]));
  // std::cout << "random seed = " << r.get_seed() << std::endl;

  // Test<CGAL::Simple_cartesian<double> >(r).run();
  Test<CGAL::Exact_predicates_inexact_constructions_kernel>(r).run();
  // Test<CGAL::Exact_predicates_exact_constructions_kernel>(r).run();

  const std::string f1 = (argc>2) ? argv[2] : CGAL::data_file_path("meshes/elephant.off");
  const std::string f2 = (argc>3) ? argv[3] : CGAL::data_file_path("meshes/sphere.off");
  // bench_on_data<CGAL::Exact_predicates_inexact_constructions_kernel>(f1,f2);

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
}
