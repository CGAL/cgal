#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/Orthogonal_standard_search.h>
#include <CGAL/Weighted_Minkowski_distance.h>
#include <CGAL/Euclidean_distance.h>
#include <CGAL/Fuzzy_iso_box_d.h>
#include <CGAL/Fuzzy_sphere_d.h>

#include <vector>
#include <iostream>

typedef CGAL::Cartesian<double>                         R;
typedef R::Point_2                                      Point_2;
typedef R::Segment_2                                    Segment_2;
typedef R::Iso_rectangle_2                              Iso_rectangle_2;
typedef R::FT                                           FT;
typedef R::Circle_2                                     Circle_2;

typedef CGAL::Creator_uniform_2<FT, Point_2>            Creator;
typedef CGAL::Plane_separator<FT>                       Separator;
typedef CGAL::Kd_tree_traits_point<Point_2>             Traits;
typedef CGAL::Euclidean_distance<Point_2>               Distance;
typedef CGAL::Orthogonal_standard_search<Traits>        Neighbour_search;
typedef CGAL::Fuzzy_iso_box_d<Point_2, Iso_rectangle_2> Fuzzy_box;
typedef CGAL::Fuzzy_sphere_d<Point_2>                   Fuzzy_circle;

typedef std::vector<Traits::Point>                      Vector;
typedef std::vector<Point_2>                            Query_vector;
