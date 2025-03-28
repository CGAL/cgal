#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Extreme_points_traits_adapter_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Convex_hull_3/predicates.h>

#include <CGAL/boost/graph/IO/polygon_mesh_io.h>

#include <vector>
#include <fstream>

// typedef CGAL::Exact_predicates_exact_constructions_kernel        K;
typedef CGAL::Exact_predicates_inexact_constructions_kernel      K;
// typedef CGAL::Simple_cartesian<double>                           K;
typedef K::Point_3                                               Point_3;
typedef K::Vector_3                                              Vector_3;
typedef CGAL::Surface_mesh<Point_3>                              Mesh;

void test_cube()
{
  std::vector<Point_3> cube;
  for(int x=0; x<2; ++x)
    for(int y=0; y<2; ++y)
      for(int z=0; z<2; ++z)
          cube.push_back(Point_3(x,y,z));

  std::vector<Point_3> inside(1, Point_3(0.25,0.25,0.20));
  std::vector<Point_3> outside(1, Point_3(-0.25,0.25,0.25));

  assert(CGAL::Convex_hull_3::do_intersect(cube, inside));
  assert(CGAL::Convex_hull_3::do_intersect(inside, cube));
  assert(!CGAL::Convex_hull_3::do_intersect(cube, outside));
  assert(!CGAL::Convex_hull_3::do_intersect(outside, cube));

  //Test intersection on vertex, edge, face
  for(double x=0.; x<=1.; x+=0.5)
    for(double y=0.; y<=1.; y+=0.5)
      for(double z=0.; z<=1.; z+=0.5){
        std::vector<Point_3> vertex(1, Point_3(x,y,z));
        assert(CGAL::Convex_hull_3::do_intersect(cube, vertex));
        assert(CGAL::Convex_hull_3::do_intersect(vertex, cube));
      }

  //Test between cubes
  std::vector<Point_3> cube_bis;
  for(int x=0; x<2; ++x)
    for(int y=0; y<2; ++y)
      for(int z=0; z<2; ++z)
          cube_bis.push_back(Point_3(x,y,z));

  auto transform=[](std::vector<Point_3> &cube, const CGAL::Aff_transformation_3<K> &t){
    for(auto &p: cube)
      p=t(p);
  };
  //Test their intersection for many translations
  for(double x=-1.5; x<=1.5; x+=0.5)
    for(double y=-1.5; y<=1.5; y+=0.5)
      for(double z=-1.5; z<=1.5; z+=0.5){
        CGAL::Aff_transformation_3<K> t(CGAL::TRANSLATION, Vector_3(x,y,z));
        transform(cube_bis, t);
        assert(CGAL::Convex_hull_3::do_intersect(cube, cube_bis)==((std::abs(x)<1.5 && std::abs(y)<1.5 && std::abs(z)<1.5)));
        transform(cube_bis, t.inverse());
      }
}

void test_degenerate()
{
  //Vertices
  std::vector<Point_3> origin(1, CGAL::ORIGIN);
  std::vector<Point_3> vertex1(1, Point_3(1,0,0));
  std::vector<Point_3> vertex2(1, Point_3(0,1,0));

  assert(CGAL::Convex_hull_3::do_intersect(origin, origin));
  assert(CGAL::Convex_hull_3::do_intersect(vertex1,vertex1));

  assert(!CGAL::Convex_hull_3::do_intersect(origin, vertex1));
  assert(!CGAL::Convex_hull_3::do_intersect(vertex1, origin));
  assert(!CGAL::Convex_hull_3::do_intersect(vertex1, vertex2));
  assert(!CGAL::Convex_hull_3::do_intersect(vertex2, vertex1));

  //Segments
  std::vector<Point_3> seg1({Point_3(0,0,0),Point_3(2,0,0)});
  std::vector<Point_3> seg2({Point_3(0,0,0),Point_3(0,2,0)});
  std::vector<Point_3> seg3({Point_3(0,2,0),Point_3(2,2,0)});
  std::vector<Point_3> seg4({Point_3(1,1,0),Point_3(1,3,0)});

  assert(CGAL::Convex_hull_3::do_intersect(origin, seg1));
  assert(CGAL::Convex_hull_3::do_intersect(seg1, origin));
  assert(CGAL::Convex_hull_3::do_intersect(vertex1, seg1));
  assert(CGAL::Convex_hull_3::do_intersect(seg1, vertex1));

  assert(!CGAL::Convex_hull_3::do_intersect(vertex2, seg1));
  assert(!CGAL::Convex_hull_3::do_intersect(seg1, vertex2));

  assert(CGAL::Convex_hull_3::do_intersect(seg1, seg1));
  assert(CGAL::Convex_hull_3::do_intersect(seg1, seg2));
  assert(CGAL::Convex_hull_3::do_intersect(seg2, seg1));
  assert(CGAL::Convex_hull_3::do_intersect(seg3, seg4));
  assert(CGAL::Convex_hull_3::do_intersect(seg4, seg3));

  assert(!CGAL::Convex_hull_3::do_intersect(seg1, seg3));
  assert(!CGAL::Convex_hull_3::do_intersect(seg3, seg1));
  assert(!CGAL::Convex_hull_3::do_intersect(seg1, seg4));
  assert(!CGAL::Convex_hull_3::do_intersect(seg4, seg1));

  //Triangle
  std::vector<Point_3> tr1({Point_3(0,0,0),Point_3(2,0,0),Point_3(0,2,0)});
  std::vector<Point_3> tr2({Point_3(0,0,0),Point_3(2,0,0),Point_3(0,0,2)});
  std::vector<Point_3> tr3({Point_3(0,0,0),Point_3(2,0,2),Point_3(0,0,2)});
  std::vector<Point_3> tr4({Point_3(1,0,0),Point_3(2,0,2),Point_3(0,0,2)});
  std::vector<Point_3> tr5({Point_3(0,0,2),Point_3(2,0,2),Point_3(0,2,2)});

  assert(CGAL::Convex_hull_3::do_intersect(tr1, tr1));
  assert(CGAL::Convex_hull_3::do_intersect(tr1, tr2));
  assert(CGAL::Convex_hull_3::do_intersect(tr2, tr1));
  assert(CGAL::Convex_hull_3::do_intersect(tr1, tr3));
  assert(CGAL::Convex_hull_3::do_intersect(tr3, tr1));
  assert(CGAL::Convex_hull_3::do_intersect(tr1, tr4));
  assert(CGAL::Convex_hull_3::do_intersect(tr4, tr1));
  assert(CGAL::Convex_hull_3::do_intersect(tr3, tr4));
  assert(CGAL::Convex_hull_3::do_intersect(tr4, tr3));

  assert(!CGAL::Convex_hull_3::do_intersect(tr1, tr5));
  assert(!CGAL::Convex_hull_3::do_intersect(tr5, tr1));
}

void test_half_sphere()
{
  std::vector<Point_3> half_sphere;
  // constexpr K::FT eps(std::pow(2,-40));
  constexpr double pi=3.14159265358979323846;
  for(double phi=25./16.; phi>0; phi-=1./4.)
    for(double theta=0; theta<2*pi; theta+=0.25)
          half_sphere.push_back(Point_3(std::sin(phi) * std::cos(theta),
                                        std::sin(phi) * std::sin(theta),
                                        std::cos(phi)));

  for(double x=-1; x<=1; x+=0.1)
    for(double y=-1; y<=1; y+=0.1){
      std::vector<Point_3> outside(1, Point_3(x,y,std::nextafter(std::cos(25./16),0)));
      std::vector<Point_3> inside(1, Point_3(x,y,std::nextafter(std::cos(25./16),1)));

      assert(CGAL::Convex_hull_3::do_intersect(half_sphere, inside)==(x*x+y*y<0.99));
      assert(!CGAL::Convex_hull_3::do_intersect(half_sphere, outside));
    }
}

int main()
{
  std::cout << std::setprecision(17);
  test_degenerate();
  test_cube();
  test_half_sphere();
  return 0;
}
