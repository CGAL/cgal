#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/Orthogonal_standard_search.h>
#include <CGAL/Weighted_Minkowski_distance.h>

#include <vector>
#include <iostream>

typedef CGAL::Cartesian<double>                         R;
typedef CGAL::Point_2<R>                                Point;
typedef CGAL::Segment_2<R>                              Segment_2;

typedef CGAL::Creator_uniform_2<double,Point>           Creator;
typedef CGAL::Plane_separator<double>                   Separator;
typedef CGAL::Kd_tree_traits_point<Separator,Point>     Traits;
typedef CGAL::Weighted_Minkowski_distance<Point, Point> Distance;
typedef CGAL::Orthogonal_standard_search<Traits, Point, Distance>
                                                        Neighbour_search;

typedef std::vector<Traits::Item>                       Vector;
typedef std::vector<Point>                              Query_vector;
