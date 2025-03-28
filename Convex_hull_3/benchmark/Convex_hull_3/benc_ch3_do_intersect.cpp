#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Extreme_points_traits_adapter_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Convex_hull_3/predicates.h>

#include <CGAL/boost/graph/IO/polygon_mesh_io.h>

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
    std::cout << t.time() << " sec" << std::endl;

  }

public:
  void run()
  {
    std::cout << "Kernel: " << typeid(K).name() << std::endl;
    Tet_tet(1000);
    Tet_mirror(1000);
    Tet_shift(1000);
    std::cout << std::endl;
  }
};

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  std::cout << "3D Distance tests" << std::endl;

  CGAL::Random rp;
  CGAL::Random r(argc==1?rp.get_seed():std::stoi(argv[1]));
  std::cout << "random seed = " << r.get_seed() << std::endl;

  // Test<CGAL::Simple_cartesian<double> >(r).run();
  // Test<CGAL::Exact_predicates_inexact_constructions_kernel>(r).run();
  Test<CGAL::Exact_predicates_exact_constructions_kernel>(r).run();

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
}
