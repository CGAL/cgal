#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/intersections.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/function_objects.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/algorithm.h>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   Kernel;
typedef Kernel::Point_3                                       Point_3;
typedef Kernel::Triangle_3                                    Triangle_3;
typedef std::vector<Triangle_3>                               Triangles;
typedef Triangles::iterator                                   Iterator;
typedef CGAL::Box_intersection_d::Box_with_handle_d<double,3,Iterator> Box;

Triangles triangles; // global vector of all triangles

// callback function that reports all truly intersecting triangles
void report_inters( const Box& a, const Box& b) {
    std::cout << "Box " << (a.handle() - triangles.begin()) << " and "
              << (b.handle() - triangles.begin()) << " intersect";
    if ( ! a.handle()->is_degenerate() && ! b.handle()->is_degenerate()
         && CGAL::do_intersect( *(a.handle()), *(b.handle()))) {
        std::cout << ", and the triangles intersect also";
    }
    std::cout << '.' << std::endl;
}

int main() {
    // Create 10 random triangles
    typedef CGAL::Random_points_in_cube_3<Point_3>           Pts;
    typedef CGAL::Creator_uniform_3< Point_3, Triangle_3>    Creator;
    typedef CGAL::Join_input_iterator_3<Pts,Pts,Pts,Creator> Triangle_gen;
    Pts    points( 1); // in centered cube [-1,1)^3
    Triangle_gen triangle_gen( points, points, points);
    CGAL::cpp11::copy_n( triangle_gen, 10, std::back_inserter(triangles));

    // Create the corresponding vector of bounding boxes
    std::vector<Box> boxes;
    for ( Iterator i = triangles.begin(); i != triangles.end(); ++i)
        boxes.push_back( Box( i->bbox(), i));

    // Run the self intersection algorithm with all defaults
    CGAL::box_self_intersection_d( boxes.begin(), boxes.end(), report_inters);
    return 0;
}
