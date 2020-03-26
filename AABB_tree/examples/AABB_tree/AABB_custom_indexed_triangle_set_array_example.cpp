#include <iostream>
#include <boost/iterator/iterator_adaptor.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>



typedef CGAL::Simple_cartesian<double> K;


// The points are stored in a flat array of doubles
// The triangles are stored in a flat array of indices
// referring to an array of coordinates: three consecutive
// coordinates represent a point, and three consecutive
// indices represent a triangle.

typedef size_t* Point_index_iterator;

// Let us now define the iterator on triangles that the tree needs:
class Triangle_iterator
    : public boost::iterator_adaptor<
    Triangle_iterator               // Derived
    , Point_index_iterator            // Base
    , boost::use_default              // Value
    , boost::forward_traversal_tag    // CategoryOrTraversal
    >
{
public:
    Triangle_iterator()
        : Triangle_iterator::iterator_adaptor_() {}

    explicit Triangle_iterator(Point_index_iterator p)
        : Triangle_iterator::iterator_adaptor_(p) {}

private:
    friend class boost::iterator_core_access;
    void increment() { this->base_reference() += 3; }
};


// The following primitive provides the conversion facilities between
// my own triangle and point types and the CGAL ones
struct My_triangle_primitive {
public:
    typedef Triangle_iterator    Id;

    // the CGAL types returned
    typedef K::Point_3    Point;
    typedef K::Triangle_3 Datum;

    // a static pointer to the vector containing the points
    // is needed to build the triangles on the fly:
    static const double* point_container;

private:
    Id m_it; // this is what the AABB tree stores internally

public:
    My_triangle_primitive() {} // default constructor needed

    // the following constructor is the one that receives the iterators from the
    // iterator range given as input to the AABB_tree
    My_triangle_primitive(Triangle_iterator a)
        : m_it(a) {}

    Id id() const { return m_it; }

    // on the fly conversion from the internal data to the CGAL types
    Datum datum() const
    {
        Point_index_iterator p_it = m_it.base();
        Point p(*(point_container + 3 * (*p_it)),
                *(point_container + 3 * (*p_it) + 1),
                *(point_container + 3 * (*p_it) + 2) );
        ++p_it;
        Point q(*(point_container + 3 * (*p_it)),
                *(point_container + 3 * (*p_it) + 1),
                *(point_container + 3 * (*p_it) + 2));
        ++p_it;
        Point r(*(point_container + 3 * (*p_it)),
                *(point_container + 3 * (*p_it) + 1),
                *(point_container + 3 * (*p_it) + 2));

        return Datum(p, q, r); // assembles triangle from three points
    }

    // one point which must be on the primitive
    Point reference_point() const
    {
      return Point(*(point_container + 3 * (*m_it)),
                   *(point_container + 3 * (*m_it) + 1),
                   *(point_container + 3 * (*m_it) + 2));
    }
};


// types
typedef CGAL::AABB_traits<K, My_triangle_primitive> My_AABB_traits;
typedef CGAL::AABB_tree<My_AABB_traits> Tree;
const double* My_triangle_primitive::point_container = 0;

int main()
{
    // generates point set
    double points[12];
    My_triangle_primitive::point_container = points;
    points[0] = 1.0; points[1] = 0.0; points[2] = 0.0;
    points[3] = 0.0; points[4] = 1.0; points[5] = 0.0;
    points[6] = 0.0; points[7] = 0.0; points[8] = 1.0;
    points[9] = 0.0; points[10] = 0.0; points[11] = 0.0;


    // generates indexed triangle set
    size_t triangles[9];
    triangles[0] = 0; triangles[1] = 1; triangles[2] = 2;
    triangles[3] = 0; triangles[4] = 1; triangles[5] = 3;
    triangles[6] = 0; triangles[7] = 3; triangles[8] = 2;

    // constructs AABB tree
    Tree tree(Triangle_iterator(triangles),
        Triangle_iterator(triangles+9));

    // counts #intersections
    K::Ray_3 ray_query(K::Point_3(0.2, 0.2, 0.2), K::Point_3(0.0, 1.0, 0.0));
    std::cout << tree.number_of_intersected_primitives(ray_query)
        << " intersections(s) with ray query" << std::endl;

    // computes closest point
    K::Point_3 point_query(2.0, 2.0, 2.0);
    K::Point_3 closest_point = tree.closest_point(point_query);
    std::cout << "closest point to " << point_query << " is: " << closest_point.x() << " " << closest_point.y() << " " << closest_point.z() << std::endl;

    return EXIT_SUCCESS;
}


