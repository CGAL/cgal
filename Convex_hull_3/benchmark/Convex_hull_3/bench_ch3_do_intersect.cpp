#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Extreme_points_traits_adapter_3.h>
#include <CGAL/convex_hull_3.h>
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

public:
   Test(CGAL::Random& r) : r(r) { }

private:
  P random_point() const
  {
    return P(FT(r.get_double(m, M)), FT(r.get_double(m, M)), FT(r.get_double(m, M)));
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

      CGAL::Convex_hull_3::predicates_impl::sphericalDisjoint(a, b, 0);
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

      CGAL::Convex_hull_3::predicates_impl::sphericalDisjoint(a, b, 0);
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

      CGAL::Convex_hull_3::predicates_impl::sphericalDisjoint(a, b, 0);
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

      CGAL::Convex_hull_3::predicates_impl::sphericalDisjoint(a, b, 0);
    }
    t.stop();
    std::cout << t.time() << " sec\n" << std::endl;

  }

public:
  void run()
  {
    std::cout << "Kernel: " << typeid(K).name() << std::endl;
    Tet_tet(1000);
    Tet_mirror(1000);
    Tet_shift(1000);
    Tet_stretched(1000);
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
  std::array<Point_3, 8> obb1, obb2;
  CGAL::oriented_bounding_box(hull1, obb1, CGAL::parameters::use_convex_hull(false));
  CGAL::oriented_bounding_box(hull2, obb2, CGAL::parameters::use_convex_hull(false));
  t.stop();
  std::cout << "Computing Obbs: " << t.time() << " sec\n" << std::endl;
  t.reset();

  t.start();
  CGAL::Convex_hull_3::do_intersect<K>(hull1, hull2);
  t.stop();
  std::cout << "Do intersect: " << t.time() << " sec\n" << std::endl;
}

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  std::cout << "3D Distance tests" << std::endl;

  CGAL::Random rp;
  // CGAL::Random r(argc==1?rp.get_seed():std::stoi(argv[1]));
  // std::cout << "random seed = " << r.get_seed() << std::endl;

  // Test<CGAL::Simple_cartesian<double> >(r).run();
  // Test<CGAL::Exact_predicates_inexact_constructions_kernel>(r).run();
  // Test<CGAL::Exact_predicates_inexact_constructions_kernel>(r).run();

  const std::string f1 = (argc>1) ? argv[1] : CGAL::data_file_path("meshes/elephant.off");
  const std::string f2 = (argc>2) ? argv[2] : CGAL::data_file_path("meshes/sphere.off");
  bench_on_data<CGAL::Exact_predicates_inexact_constructions_kernel>(f1,f2);

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
}
