// Author(s): Camille Wormser, Pierre Alliez
// Example of an AABB tree used with indexed triangle set

#include <iostream>
#include <list>
#include <boost/iterator.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>



typedef CGAL::Simple_cartesian<double> K;


// My own point type:
struct My_point {
    double x;
    double y;
    double z;

    My_point (double _x, double _y, double _z)
        : x(_x), y(_y), z(_z) {}
};

// The triangles are stored in a flat array of indices 
// referring to an array of points: three consecutive
// indices represent a triangle.
typedef std::vector<size_t>::const_iterator Point_index_iterator;

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
    static const std::vector<My_point>* point_container;

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
        const My_point& mp = (*point_container)[*p_it];
        Point p(mp.x, mp.y, mp.z);
        ++p_it;
        const My_point& mq = (*point_container)[*p_it];
        Point q(mq.x, mq.y, mq.z);
        ++p_it;
        const My_point& mr = (*point_container)[*p_it];
        Point r(mr.x, mr.y, mr.z);

        return Datum(p, q, r); // assembles triangle from three points
    }

    // one point which must be on the primitive
    Point reference_point() const
    { 
        const My_point& mp = (*point_container)[*m_it];
        return Point(mp.x, mp.y, mp.z);
    }
};


// types
typedef CGAL::AABB_traits<K, My_triangle_primitive> My_AABB_traits;
typedef CGAL::AABB_tree<My_AABB_traits> Tree;
const std::vector<My_point>* My_triangle_primitive::point_container = 0;

int main()
{
    // generates point set
    My_point a(1.0, 0.0, 0.0);
    My_point b(0.0, 1.0, 0.0);
    My_point c(0.0, 0.0, 1.0);
    My_point d(0.0, 0.0, 0.0);

    std::vector<My_point> points;
    My_triangle_primitive::point_container = &points;
    points.push_back(a);
    points.push_back(b);
    points.push_back(c);
    points.push_back(d); 

    // generates indexed triangle set
    std::vector<size_t> triangles;
    triangles.push_back(0); triangles.push_back(1); triangles.push_back(2);  
    triangles.push_back(0); triangles.push_back(1); triangles.push_back(3);  
    triangles.push_back(0); triangles.push_back(3); triangles.push_back(2);  

    // constructs AABB tree
    Tree tree(Triangle_iterator(triangles.begin()), 
        Triangle_iterator(triangles.end()));

    // counts #intersections
    K::Ray_3 ray_query(K::Point_3(0.2, 0.2, 0.2), K::Point_3(0.0, 1.0, 0.0));
    std::cout << tree.number_of_intersected_primitives(ray_query)
        << " intersections(s) with ray query" << std::endl;

    // computes closest point
    K::Point_3 point_query(2.0, 2.0, 2.0);
    K::Point_3 closest_point = tree.closest_point(point_query);

    std::cerr << "closest point is: " << closest_point << std::endl;
    return EXIT_SUCCESS;
}
