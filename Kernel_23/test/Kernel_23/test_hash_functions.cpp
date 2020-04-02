#include <unordered_set>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Bbox_2.h>

typedef CGAL::Simple_cartesian<double> SC;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Epick;

template <typename Object>
void test (const Object& obj)
{
  std::unordered_set<Object> unordered_set;
  unordered_set.insert (obj);
  unordered_set.insert (obj);
  assert (unordered_set.size() == 1);
}

int main()
{
  test (SC::Aff_transformation_2 (1, 2, 3, 4));
  test (Epick::Aff_transformation_2 (1, 2, 3, 4));

  test (CGAL::Bbox_2 ());

  test (SC::Circle_2 (CGAL::ORIGIN, 1));
  test (Epick::Circle_2 (CGAL::ORIGIN, 1));

  test (SC::Iso_rectangle_2 (1, 2, 3, 4));
  test (Epick::Iso_rectangle_2 (1, 2, 3, 4));

  test (SC::Point_2 (1, 2));
  test (Epick::Point_2 (1, 2));

  test (SC::Segment_2 (SC::Point_2 (1, 2), SC::Point_2 (3, 4)));
  test (Epick::Segment_2 (Epick::Point_2 (1, 2), Epick::Point_2 (3, 4)));

  test (SC::Vector_2 (1, 2));
  test (Epick::Vector_2 (1, 2));

  test (SC::Weighted_point_2 (SC::Point_2 (1, 2), 3));
  test (Epick::Weighted_point_2 (Epick::Point_2 (1, 2), 3));

  test (SC::Aff_transformation_3 (1, 2, 3, 4, 5, 6, 7, 8, 9));
  test (Epick::Aff_transformation_3 (1, 2, 3, 4, 5, 6, 7, 8, 9));

  test (CGAL::Bbox_3 ());

  test (SC::Iso_cuboid_3 (1, 2, 3, 4, 5, 6));
  test (Epick::Iso_cuboid_3 (1, 2, 3, 4, 5, 6));

  test (SC::Point_3 (1, 2, 3));
  test (Epick::Point_3 (1, 2, 3));

  test (SC::Segment_3 (SC::Point_3 (1, 2, 3), SC::Point_3 (4, 5, 6)));
  test (Epick::Segment_3 (Epick::Point_3 (1, 2, 3), Epick::Point_3 (4, 5, 6)));

  test (SC::Sphere_3 (CGAL::ORIGIN, 1));
  test (Epick::Sphere_3 (CGAL::ORIGIN, 1));

  test (SC::Vector_3 (1, 2, 3));
  test (Epick::Vector_3 (1, 2, 3));

  test (SC::Weighted_point_3 (SC::Point_3 (1, 2, 3), 4));
  test (Epick::Weighted_point_3 (Epick::Point_3 (1, 2, 3), 4));

  return 0;
}


