#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_set_2.h>
#include <list>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;


int main()
{
  CGAL::Polygon_2<K> ob;
  ob.push_back(CGAL::Point_2<K>(1, 1));
  ob.push_back(CGAL::Point_2<K>(1, 0));
  ob.push_back(CGAL::Point_2<K>(6, 0));
  ob.push_back(CGAL::Point_2<K>(6, 7));
  ob.push_back(CGAL::Point_2<K>(0, 7));
  ob.push_back(CGAL::Point_2<K>(0, 1));

  CGAL::Polygon_2<K> h;
  h.push_back(CGAL::Point_2<K>(2, 1));
  h.push_back(CGAL::Point_2<K>(2, 2));
  h.push_back(CGAL::Point_2<K>(3, 2));
  h.push_back(CGAL::Point_2<K>(3, 3));
  h.push_back(CGAL::Point_2<K>(2, 3));
  h.push_back(CGAL::Point_2<K>(2, 4));
  h.push_back(CGAL::Point_2<K>(3, 4));
  h.push_back(CGAL::Point_2<K>(3, 5));
  h.push_back(CGAL::Point_2<K>(4, 5));
  h.push_back(CGAL::Point_2<K>(4, 1));

  CGAL::Polygon_with_holes_2<K> ob_with_holes(ob);
  ob_with_holes.add_hole(h);
  CGAL::Polygon_set_2<K> inter(ob_with_holes);

  CGAL::Polygon_2<K> new_poly;
  new_poly.push_back(CGAL::Point_2<K>(1, 1));
  new_poly.push_back(CGAL::Point_2<K>(2, 1));
  new_poly.push_back(CGAL::Point_2<K>(2, 2));
  new_poly.push_back(CGAL::Point_2<K>(2, 3));
  new_poly.push_back(CGAL::Point_2<K>(2, 4));
  new_poly.push_back(CGAL::Point_2<K>(2, 5));
  new_poly.push_back(CGAL::Point_2<K>(3, 5));
  new_poly.push_back(CGAL::Point_2<K>(4, 5));
  new_poly.push_back(CGAL::Point_2<K>(4, 6));
  new_poly.push_back(CGAL::Point_2<K>(1, 6));

  inter.difference(new_poly);
}