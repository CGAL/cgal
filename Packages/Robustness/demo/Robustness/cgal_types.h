#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>

#include <CGAL/segment_intersection_points_2.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/copy_n.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Cartesian_converter.h> 
#include <CGAL/kernel_to_kernel.h>
#include <CGAL/Timer.h>


typedef CGAL::Quotient<CGAL::MP_Float>  exact_NT;
typedef CGAL::Cartesian<float>          CartesianFloat;
typedef CGAL::Cartesian<double>         CartesianDouble;
typedef CGAL::Homogeneous<float>        HomogeneousFloat;
typedef CGAL::Homogeneous<double>       HomogeneousDouble;

typedef CGAL::Cartesian<double>         C_double;
typedef C_double::Point_2               double_Point;
typedef C_double::Segment_2             double_Segment;
typedef CGAL::Cartesian<exact_NT>       C_real;
typedef C_real::Point_2                 real_Point;
typedef C_real::Segment_2               real_Segment;
typedef CGAL::Creator_uniform_2<double, double_Point>
                                        Point_creator;
typedef CGAL::Random_points_in_square_2<double_Point, Point_creator>
                                        Source;
typedef CGAL::Creator_uniform_2<double_Point,  double_Segment>
                                        Segment_creator;
typedef CGAL::Join_input_iterator_2<Source, Source, Segment_creator>
                                        Segment_iterator;
