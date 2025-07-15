#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/Extreme_points_traits_adapter_3.h>

#include <CGAL/Convex_hull_3/predicates.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Convex_hull_hierarchy.h>

#include <boost/property_map/vector_property_map.hpp>

#include <vector>
#include <fstream>

typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt        EPECK_with_sqrt;
typedef CGAL::Exact_predicates_exact_constructions_kernel        EPECK;
typedef CGAL::Exact_predicates_inexact_constructions_kernel      EPICK;
typedef CGAL::Simple_cartesian<double>                           DOUBLE;

const double PI=3.14159265358979323846;

template<typename K>
struct Sphere{
  typedef typename K::FT FT;
  typedef typename K::Point_3 P;

  //do_intersect deduce Kernel from iterator type of the object, this an astuce for kernel deduction
  typedef typename std::vector<P>::iterator iterator;

  FT r;
  P c;
  Sphere(P c_, FT r_):r(r_), c(c_){}

  template<class Vector_3, class Converter>
  P extreme_point(const Vector_3 dir, const Converter& converter) const{
    return converter(c)+converter(r)*dir/CGAL::sqrt(dir.squared_length());
  }
};
template<typename K, class Vector_3, class Converter>
typename K::Point_3 extreme_point(const Sphere<K> &sp, const Vector_3 &dir, const Converter& converter){
  return sp.extreme_point(dir, converter);
}

template<typename K, typename Mesh>
struct Test{
  typedef typename K::Point_3                                      P;
  typedef typename K::Vector_3                                     V;
  typedef typename K::FT                                           FT;

  typedef boost::vector_property_map<P> PMap;

  void test(std::vector<P> &vec_a, const std::vector<P> &vec_b, bool result){
    assert(CGAL::Convex_hull_3::do_intersect(vec_a, vec_b)==result);
    assert(CGAL::Convex_hull_3::do_intersect(vec_b, vec_a)==result);

    Mesh sm_a, sm_b;
    CGAL::convex_hull_3(vec_a.begin(), vec_a.end(), sm_a);
    CGAL::convex_hull_3(vec_b.begin(), vec_b.end(), sm_b);
    assert(CGAL::Convex_hull_3::do_intersect(sm_a, sm_b)==result);
    assert(CGAL::Convex_hull_3::do_intersect(sm_b, sm_a)==result);

    if constexpr(std::is_same_v<Mesh, CGAL::Surface_mesh<P> >){
      CGAL::Convex_hull_hierarchy<Mesh> hsm_a(vec_a.begin(),vec_a.end()), hsm_b(vec_b.begin(),vec_b.end());
      assert(CGAL::Convex_hull_3::do_intersect(hsm_a, hsm_b)==result);
      assert(CGAL::Convex_hull_3::do_intersect(hsm_b, hsm_a)==result);

      CGAL::Convex_hull_hierarchy<Mesh> hsm_a_2(sm_a), hsm_b_2(sm_b);
      assert(CGAL::Convex_hull_3::do_intersect(hsm_a_2, hsm_b_2)==result);
      assert(CGAL::Convex_hull_3::do_intersect(hsm_b_2, hsm_a_2)==result);
    }

    //Test with Point map
    PMap v_a, v_b; //Vectors de Point_3
    for(size_t i=0; i<vec_a.size(); ++i)
      v_a[i]=vec_a[i];
    for(size_t i=0; i<vec_b.size(); ++i)
      v_b[i]=vec_b[i];
    std::vector< size_t > pm_a(vec_a.size(),0), pm_b(vec_b.size(), 0); //Vector d'indices
    std::iota(pm_a.begin(), pm_a.end(), 0);
    std::iota(pm_b.begin(), pm_b.end(), 0);
    assert(CGAL::Convex_hull_3::do_intersect(pm_a, pm_b, CGAL::parameters::point_map(v_a), CGAL::parameters::point_map(v_b))==result);
    assert(CGAL::Convex_hull_3::do_intersect(pm_b, pm_a, CGAL::parameters::point_map(v_b), CGAL::parameters::point_map(v_a))==result);

    CGAL::Surface_mesh<size_t> sm_a_pm, sm_b_pm;
    CGAL::convex_hull_3(pm_a.begin(), pm_a.end(), sm_a_pm, CGAL::make_extreme_points_traits_adapter(v_a));
    CGAL::convex_hull_3(pm_b.begin(), pm_b.end(), sm_b_pm, CGAL::make_extreme_points_traits_adapter(v_b));
    assert(CGAL::Convex_hull_3::do_intersect(sm_a_pm, sm_b_pm, CGAL::parameters::vertex_point_map(make_compose_property_map(sm_a_pm.points(), v_a)),
                                                             CGAL::parameters::vertex_point_map(make_compose_property_map(sm_b_pm.points(), v_b)))==result);
    assert(CGAL::Convex_hull_3::do_intersect(sm_b_pm, sm_a_pm, CGAL::parameters::vertex_point_map(make_compose_property_map(sm_b_pm.points(), v_b)),
                                                             CGAL::parameters::vertex_point_map(make_compose_property_map(sm_a_pm.points(), v_a)))==result);

    if constexpr(std::is_same_v<Mesh, CGAL::Surface_mesh<P> >){
      CGAL::Surface_mesh< typename Mesh::Vertex_index > sma_pm2, smb_pm2;
      CGAL::convex_hull_3(vertices(sm_a).begin(), vertices(sm_a).end(), sma_pm2, CGAL::make_extreme_points_traits_adapter(sm_a.points()));
      CGAL::convex_hull_3(vertices(sm_b).begin(), vertices(sm_b).end(), smb_pm2, CGAL::make_extreme_points_traits_adapter(sm_b.points()));
      assert(CGAL::Convex_hull_3::do_intersect(sma_pm2, smb_pm2, CGAL::parameters::vertex_point_map(make_compose_property_map(sma_pm2.points(), sm_a.points())),
                                                                CGAL::parameters::vertex_point_map(make_compose_property_map(smb_pm2.points(), sm_b.points())))==result);
      assert(CGAL::Convex_hull_3::do_intersect(smb_pm2, sma_pm2, CGAL::parameters::vertex_point_map(make_compose_property_map(smb_pm2.points(), sm_b.points())),
                                                                CGAL::parameters::vertex_point_map(make_compose_property_map(sma_pm2.points(), sm_a.points())))==result);


      // CGAL::Convex_hull_hierarchy<size_t> hsma_pm(sma_pm, CGAL::make_extreme_points_traits_adapter(va));
      // CGAL::Convex_hull_with_hierarchy<size_t> hsmb_pm(smb_pm, CGAL::make_extreme_points_traits_adapter(vb));
      // assert(CGAL::Convex_hull_3::do_intersect(hsma_pm, hsmb_pm, CGAL::Convex_hull_3::make_do_intersect_traits_with_point_maps(va, vb))==result);
      // assert(CGAL::Convex_hull_3::do_intersect(hsmb_pm, hsma_pm, CGAL::Convex_hull_3::make_do_intersect_traits_with_point_maps(vb, va))==result);

      // CGAL::Convex_hull_hierarchy<CGAL::Surface_mesh< typename Mesh::Vertex_index > > hsma_pm2(sm_a, CGAL::make_extreme_points_traits_adapter(sm_a.points()));
      // CGAL::Convex_hull_hierarchy<CGAL::Surface_mesh< typename Mesh::Vertex_index > > hsmb_pm2(sm_b, CGAL::make_extreme_points_traits_adapter(sm_b.points()));
      // assert(CGAL::Convex_hull_3::do_intersect(hsma_pm2, hsmb_pm2, CGAL::parameters::vertex_point_map(sm_a.points()), CGAL::parameters::vertex_point_map(sm_b.points()))==result);
      // assert(CGAL::Convex_hull_3::do_intersect(hsmb_pm2, hsma_pm2, CGAL::parameters::vertex_point_map(sm_b.points()), CGAL::parameters::vertex_point_map(sm_a.points()))==result);
    }
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


    // constexpr double pi=3.14159265358979323846;
    // for(double phi=25./16.; phi>0; phi-=1./4.)
    //   for(double theta=0; theta<2*pi; theta+=0.25)
    //         half_sphere.push_back(P(std::sin(phi) * std::cos(theta),
    //                                 std::sin(phi) * std::sin(theta),
    //                                 std::cos(phi)));

    for(double x=-0.5; x<=0.5; x+=0.1)
      for(double y=-0.5; y<=0.5; y+=0.1){
        std::vector<P> outside(1, P(x,y,std::nextafter(std::cos(25./16),0)));
        std::vector<P> inside(1, P(x,y,std::nextafter(std::cos(25./16),1)));

        test(half_sphere, inside, true);
        test(half_sphere, outside, false);
      }
  }

  void test_random_tetrahedra(int N, CGAL::Random &r)
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
      test(a, b, CGAL::do_intersect(Tet(p0, p1, p2, p3), Tet(q0, q1, q2, q3)));
    }
  }

  void test_random_sphere(int N, int M, CGAL::Random &r)
  {
    auto random_sphere_point=[&](){
      double a=r.get_double(0,2*PI);
      double b=r.get_double(-PI/2,PI/2);
      return P( std::cos(a)*std::cos(b), std::sin(a)*std::cos(b), std::sin(b));
    };
    for(int i=0; i<M; ++i)
    {
      std::vector<P> sp_a;
      std::vector<P> sp_b;
      V vec1( P(0,0,0), random_sphere_point());
      V vec2( P(0,0,0), random_sphere_point());
      vec1*=0.9;
      vec2*=0.9;
      for(int i=0; i<N; ++i){
        sp_a.push_back(random_sphere_point()+vec1);
        sp_b.push_back(random_sphere_point()+vec2);
      }
      test(sp_a, sp_b, true);
    }
  }

  void test_implicit_function(CGAL::Random &r, int N){
    std::vector< Sphere<K> > spheres;
    for(int i=0; i<N; ++i)
      spheres.emplace_back(P(r.get_double(-2, 2),0,0),1);

    for(int i=0; i<N; ++i)
      for(int j=0; j<N; ++j)
        assert(CGAL::Convex_hull_3::do_intersect(spheres[i], spheres[j])==(CGAL::abs(spheres[i].c.x()-spheres[j].c.x())<=2));
  }

  void full_test(CGAL::Random &r){
    // test_degenerate();
    test_cube();
    test_half_sphere();
    test_random_tetrahedra(1000, r);
    test_random_sphere(5000,20, r);
  }

};

int main(int argc, char** argv)
{
  CGAL::Random r(argc==1?CGAL::get_default_random():std::stoi(argv[1]));
  std::cout << "random seed = " << r.get_seed() << std::endl;

  std::cout << std::setprecision(17);
  // Test<DOUBLE>().full_test(r);
  Test<EPICK, CGAL::Surface_mesh<EPICK::Point_3> >().full_test(r);
  Test<EPICK, CGAL::Polyhedron_3<EPICK> >().full_test(r);
  Test<EPECK, CGAL::Surface_mesh<EPECK::Point_3> >().full_test(r);
  // Test<EPECK_with_sqrt>().test_implicit_function(r,10);
  return 0;
}
