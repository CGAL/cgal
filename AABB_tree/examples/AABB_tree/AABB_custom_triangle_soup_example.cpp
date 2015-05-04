// Author(s) : Camille Wormser, Pierre Alliez
// Example of an AABB tree used with a simple list of 
// triangles (a triangle soup) stored into an array of points.

#include <iostream>
#include <vector>
#include <boost/iterator.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>



typedef CGAL::Simple_cartesian<double> K;


// My own point type
struct My_point {
    double x;
    double y;
    double z;

    My_point (double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
};

// The triangles are stored in a flat vector of points (a triangle soup): 
// three consecutive points represent a triangle
typedef std::vector<My_point>::const_iterator Point_iterator;

// defines the iterator over triangles needed by the tree:
class Triangle_iterator
    : public boost::iterator_adaptor<
    Triangle_iterator               // Derived
    , Point_iterator                  // Base
    , boost::use_default              // Value
    , boost::forward_traversal_tag    // CategoryOrTraversal
    >
{
public:
    Triangle_iterator()
        : Triangle_iterator::iterator_adaptor_() {}

    explicit Triangle_iterator(Point_iterator p)
        : Triangle_iterator::iterator_adaptor_(p) {}

private:
    friend class boost::iterator_core_access;
    void increment() { this->base_reference() += 3; }
};


// The following primitive provides the conversion facilities between
// my own triangle and point types and the CGAL ones
struct My_triangle_primitive {
public:
    typedef Triangle_iterator Id;

    // the CGAL types returned
    typedef K::Point_3    Point;
    typedef K::Triangle_3 Datum;
private:
    Id m_it; // this is what the AABB tree will store internally

public:
    My_triangle_primitive() {} // default constructor needed

    // the following constructor is the one that receives the iterators from the 
    // iterator range given as input to the AABB_tree
    My_triangle_primitive(Triangle_iterator a)
        : m_it(a) {}

    Id id() const { return m_it; }

    // on the fly conversion from the internal data
    // to the CGAL types
    Datum datum() const
    { 
        Point_iterator p_it = m_it.base();
        Point p(p_it->x, p_it->y, p_it->z);
        ++p_it;
        Point q(p_it->x, p_it->y, p_it->z);
        ++p_it;
        Point r(p_it->x, p_it->y, p_it->z);

        return Datum(p, q, r); // assembles a triangle from three points
    }

    // returns one point which must be on the primitive
    Point reference_point() const
    { 
        return Point(m_it->x, m_it->y, m_it->z);
    }
};

// types
typedef CGAL::AABB_traits<K, My_triangle_primitive> My_AABB_traits;
typedef CGAL::AABB_tree<My_AABB_traits> Tree;

int main()
{
    // generates triangle soup
    My_point a(1.0, 0.0, 0.0);
    My_point b(0.0, 1.0, 0.0);
    My_point c(0.0, 0.0, 1.0);
    My_point d(0.0, 0.0, 0.0);

    std::vector<My_point> triangles;
    triangles.push_back(a); triangles.push_back(b); triangles.push_back(c);  
    triangles.push_back(a); triangles.push_back(b); triangles.push_back(d);  
    triangles.push_back(a); triangles.push_back(d); triangles.push_back(c);  

    // constructs AABB tree
    Tree tree(Triangle_iterator(triangles.begin()), 
        Triangle_iterator(triangles.end()));

    // counts #intersections
    K::Ray_3 ray_query(K::Point_3(1.0, 0.0, 0.0), K::Point_3(0.0, 1.0, 0.0));
    std::cout << tree.number_of_intersected_primitives(ray_query)
        << " intersections(s) with ray query" << std::endl;

    // computes closest point
    K::Point_3 point_query(2.0, 2.0, 2.0);
    K::Point_3 closest_point = tree.closest_point(point_query);

    std::cerr << "closest point is: " << closest_point << std::endl;
    return EXIT_SUCCESS;
}
