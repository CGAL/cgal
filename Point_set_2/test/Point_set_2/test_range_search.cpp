#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <cassert>
#include <CGAL/Point_set_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <list>


typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef K::Point_2         Point_2;
typedef K::Triangle_2      Triangle_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;
typedef K::Circle_2        Circle_2;
typedef CGAL::Point_set_2<K>::Edge    Edge;
typedef CGAL::Point_set_2<K>::Edge_iterator  Edge_iterator;
typedef CGAL::Point_set_2<K>::Vertex_handle  Vertex_handle;
typedef CGAL::Point_set_2<K>::Vertex  Vertex;

CGAL::Point_set_2<K> PSet;


template <typename Shape>
struct Has_not_on_unbounded_side {
  Has_not_on_unbounded_side(const Shape& s)
   : shape(s)
  {}

  bool operator()(const Point_2& p) const
  {
    return ! shape.has_on_unbounded_side(p);
  }

private:
  Shape shape;
};


int main()
{
  CGAL::Random rnd = CGAL::Random(7812);
  CGAL::Random_points_in_disc_2<Point_2> rpg(100.0, rnd);
  std::list<Point_2> points;
  std::list<Vertex_handle> LV;

  std::copy_n(rpg, 1000, std::back_inserter(points));
  PSet.insert(points.begin(), points.end());

  Point_2 p(10, 10), q(50, 10), r(50, 50), s(10, 50);

  std::cout << "range search for circle" << std::endl;
  Circle_2 circle(p, q);
  PSet.range_search(circle, back_inserter(LV));

  std::list<Point_2>::difference_type m =
           std::count_if(points.begin(),
                         points.end(),
                         Has_not_on_unbounded_side<Circle_2>(circle));
  assert((std::size_t) m == LV.size());


  std::cout << "range search for triangle" << std::endl;

  LV.clear();
  PSet.range_search(p ,q, r,back_inserter(LV));

  m =  std::count_if(points.begin(),
                     points.end(),
                     Has_not_on_unbounded_side<Triangle_2>(Triangle_2(p,q,r)));
  assert((std::size_t) m == LV.size());

  std::cout << "range search for iso rectangle" << std::endl;
  LV.clear();
  PSet.range_search(p ,q, r ,s, back_inserter(LV));

  m =  std::count_if(points.begin(),
                     points.end(),
                     Has_not_on_unbounded_side<Iso_rectangle_2>(Iso_rectangle_2(p,r)));
  assert((std::size_t) m == LV.size());
  return 0;
}
