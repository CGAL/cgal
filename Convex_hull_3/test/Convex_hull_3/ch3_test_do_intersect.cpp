#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/Extreme_points_traits_adapter_3.h>

#include <CGAL/Convex_hull_3/predicates.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/convex_hull_with_hierarchy.h>

#include <vector>
#include <fstream>

typedef CGAL::Exact_predicates_exact_constructions_kernel        EPECK;
typedef CGAL::Exact_predicates_inexact_constructions_kernel      EPICK;
typedef CGAL::Simple_cartesian<double>                           DOUBLE;

template<typename K>
struct Test{
  typedef typename K::Point_3                                      P;
  typedef typename K::Vector_3                                     V;
  typedef CGAL::Surface_mesh<P>                                 Mesh;

  void test(std::vector<P> &a, std::vector<P> &b, bool result){
    assert(CGAL::Convex_hull_3::do_intersect<K>(a, b)==result);
    assert(CGAL::Convex_hull_3::do_intersect<K>(b, a)==result);
    CGAL::Surface_mesh<P> sma, smb;
    CGAL::convex_hull_3(a.begin(), a.end(), sma);
    CGAL::convex_hull_3(b.begin(), b.end(), smb);
    assert(CGAL::Convex_hull_3::do_intersect<K>(sma, smb)==result);
    assert(CGAL::Convex_hull_3::do_intersect<K>(smb, sma)==result);
    CGAL::Convex_hull_with_hierarchy<P> hsma(sma), hsmb(smb);
    assert(CGAL::Convex_hull_3::do_intersect<K>(hsma, hsmb)==result);
    assert(CGAL::Convex_hull_3::do_intersect<K>(hsmb, hsma)==result);
  }

  void test_cube()
  {
    std::vector<P> cube;
    for(int x=0; x<2; ++x)
      for(int y=0; y<2; ++y)
        for(int z=0; z<2; ++z)
            cube.push_back(P(x,y,z));

    std::vector<P> inside(1, P(0.25,0.25,0.20));
    std::vector<P> outside(1, P(-0.25,0.25,0.25));

    test(cube, inside, true);
    test(cube, outside, false);

    //Test intersection on vertex, edge, face
    for(double x=0.; x<=1.; x+=0.5)
      for(double y=0.; y<=1.; y+=0.5)
        for(double z=0.; z<=1.; z+=0.5){
          std::vector<P> vertex(1, P(x,y,z));
          test(cube, vertex, true);
        }

    //Test between cubes
    std::vector<P> cube_bis;
    for(int x=0; x<2; ++x)
      for(int y=0; y<2; ++y)
        for(int z=0; z<2; ++z)
            cube_bis.push_back(P(x,y,z));

    auto transform=[](std::vector<P> &cube, const CGAL::Aff_transformation_3<K> &t){
      for(auto &p: cube)
        p=t(p);
    };
    //Test their intersection for many translations
    for(double x=-1.5; x<=1.5; x+=0.5)
      for(double y=-1.5; y<=1.5; y+=0.5)
        for(double z=-1.5; z<=1.5; z+=0.5){
          CGAL::Aff_transformation_3<K> t(CGAL::TRANSLATION, V(x,y,z));
          transform(cube_bis, t);
          test(cube, cube_bis, std::abs(x)<1.5 && std::abs(y)<1.5 && std::abs(z)<1.5);
          transform(cube_bis, t.inverse());
        }
  }

  void test_degenerate()
  {
    //Vertices
    std::vector<P> origin(1, CGAL::ORIGIN);
    std::vector<P> vertex1(1, P(1,0,0));
    std::vector<P> vertex2(1, P(0,1,0));

    test(origin, origin, true);
    test(vertex1,vertex1, true);

    test(origin, vertex1, false);
    test(vertex1, vertex2, false);

    //Segments
    std::vector<P> seg1({P(0,0,0),P(2,0,0)});
    std::vector<P> seg2({P(0,0,0),P(0,2,0)});
    std::vector<P> seg3({P(0,2,0),P(2,2,0)});
    std::vector<P> seg4({P(1,1,0),P(1,3,0)});

    test(origin, seg1, true);
    test(vertex1, seg1, true);

    test(vertex2, seg1, false);

    test(seg1, seg1, true);
    test(seg1, seg2, true);
    test(seg3, seg4, true);

    test(seg1, seg3, false);
    test(seg1, seg4, false);

    //Triangle
    std::vector<P> tr1({P(0,0,0),P(2,0,0),P(0,2,0)});
    std::vector<P> tr2({P(0,0,0),P(2,0,0),P(0,0,2)});
    std::vector<P> tr3({P(0,0,0),P(2,0,2),P(0,0,2)});
    std::vector<P> tr4({P(1,0,0),P(2,0,2),P(0,0,2)});
    std::vector<P> tr5({P(0,0,2),P(2,0,2),P(0,2,2)});

    test(tr1, tr1, true);
    test(tr1, tr2, true);
    test(tr1, tr3, true);
    test(tr1, tr4, true);
    test(tr3, tr4, true);

    test(tr1, tr5, false);
    test(tr5, tr1, false);
  }

  void test_half_sphere()
  {
    std::vector<P> half_sphere;
    // constexpr K::FT eps(std::pow(2,-40));
    constexpr double pi=3.14159265358979323846;
    for(double phi=25./16.; phi>0; phi-=1./4.)
      for(double theta=0; theta<2*pi; theta+=0.25)
            half_sphere.push_back(P(std::sin(phi) * std::cos(theta),
                                    std::sin(phi) * std::sin(theta),
                                    std::cos(phi)));

    for(double x=-1; x<=1; x+=0.1)
      for(double y=-1; y<=1; y+=0.1){
        std::vector<P> outside(1, P(x,y,std::nextafter(std::cos(25./16),0)));
        std::vector<P> inside(1, P(x,y,std::nextafter(std::cos(25./16),1)));

        test(half_sphere, inside, x*x+y*y<0.99);
        test(half_sphere, outside, false);
      }
  }

  void test_random_tetrahedrons(int N, CGAL::Random &r)
  {
    using Tet = typename K::Tetrahedron_3;

    auto random_point=[&](){
      return P(r.get_double(0, 1), r.get_double(0, 1), r.get_double(0, 1));
    };
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

      std::vector<P> a({p1,p2,p3,p0});
      std::vector<P> b({q1,q2,q3,q0});
      test(a, b, CGAL::do_intersect<K>(Tet(p0, p1, p2, p3), Tet(q0, q1, q2, q3)));
    }
  }

  void full_test(CGAL::Random &r){
    test_degenerate();
    test_cube();
    test_half_sphere();
    test_random_tetrahedrons(1000, r);
  }

};

int main(int argc, char** argv)
{
  CGAL::Random rp;
  CGAL::Random r(argc==1?rp.get_seed():std::stoi(argv[1]));
  std::cout << "random seed = " << r.get_seed() << std::endl;

  std::cout << std::setprecision(17);
  // Test<DOUBLE>().full_test(r);
  Test<EPICK>().full_test(r);
  Test<EPECK>().full_test(r);
  return 0;
}
