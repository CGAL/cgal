
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/property_map.h>

typedef CGAL::Simple_cartesian<double> K1;
typedef CGAL::Simple_cartesian<CGAL::Quotient<CGAL::MP_Float> >  K2;
typedef K1::Point_3 Point_3;

template <typename Mesh>
void
test()
{
  Mesh m;
  CGAL::make_triangle(Point_3(2,0,0),Point_3(1,0,0),Point_3(1,1,0),m);

  typedef typename boost::property_map<Mesh, CGAL::vertex_point_t >::type VPMap;
  VPMap vmap = get(CGAL::vertex_point, m);

  CGAL::Cartesian_converter_property_map<K2::Point_3, VPMap> kcmap =CGAL::make_cartesian_converter_property_map<K2::Point_3>(vmap);
  CGAL_assertion(get(kcmap, *vertices(m).begin()) == CGAL::Point_3<K2>(2,0,0));
  put(kcmap, *vertices(m).begin(), CGAL::Point_3<K2>(0,2,3));
  CGAL_assertion(get(kcmap, *vertices(m).begin()) == CGAL::Point_3<K2>(0,2,3));

}

int main()
{

  typedef CGAL::Surface_mesh<Point_3> SM;
  test<SM>();
  return 0;
}

