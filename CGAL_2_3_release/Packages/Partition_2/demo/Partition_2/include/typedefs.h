#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Aff_transformation_2.h>
#include <list>


typedef CGAL::Cartesian<double>                           K;
typedef K::Vector_2                                       Vector_2;
typedef CGAL::Polygon_traits_2<K>                         Traits;
typedef Traits::Point_2                                   Point_2;
typedef std::list<Point_2>                                Container;
typedef CGAL::Polygon_2<Traits, Container>                Polygon_2;
typedef CGAL::Creator_uniform_2<int, Point_2>             Creator;
typedef CGAL::Random_points_in_square_2<Point_2, Creator> Point_generator;
typedef CGAL::Aff_transformation_2<K>                     Transformation_2;



#endif // TYPEDEFS_H
