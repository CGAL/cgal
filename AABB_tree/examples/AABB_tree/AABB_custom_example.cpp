// Author(s): Camille Wormser, Pierre Alliez
// An example of an AABB tree constructed with custom point and triangle types.

#include <iostream>
#include <list>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>



typedef CGAL::Simple_cartesian<double> K;


// custom point type
struct My_point {
    double m_x;
    double m_y;
    double m_z;

    My_point(const double x,
        const double y,
        const double z)
        : m_x(x), m_y(y), m_z(z) {}
};

// custom triangle type with
// three pointers to points
struct My_triangle {
    My_point *m_pa;
    My_point *m_pb;
    My_point *m_pc;

    My_triangle(My_point *pa,
        My_point *pb,
        My_point *pc)
        : m_pa(pa), m_pb(pb), m_pc(pc) {}
};

// the custom triangles are stored into a vector
typedef std::vector<My_triangle>::const_iterator Iterator;

// The following primitive provides the conversion facilities between
// the custom triangle and point types and the CGAL ones
struct My_triangle_primitive {
public:

    // this is the type of data that the queries returns. For this example
    // we imagine that, for some reasons, we do not want to store the iterators
    // of the vector, but raw pointers. This is to show that the Id type
    // does not have to be the same as the one of the input parameter of the
    // constructor.
    typedef const My_triangle* Id;

    // CGAL types returned
    typedef K::Point_3    Point; // CGAL 3D point type
    typedef K::Triangle_3 Datum; // CGAL 3D triangle type

private:
    Id m_pt; // this is what the AABB tree stores internally

public:
    My_triangle_primitive() {} // default constructor needed

    // the following constructor is the one that receives the iterators from the
    // iterator range given as input to the AABB_tree
    My_triangle_primitive(Iterator it)
        : m_pt(&(*it)) {}

    const Id& id() const { return m_pt; }

    // utility function to convert a custom
    // point type to CGAL point type.
    Point convert(const My_point *p) const
    {
        return Point(p->m_x,p->m_y,p->m_z);
    }

    // on the fly conversion from the internal data to the CGAL types
    Datum datum() const
    {
        return Datum(convert(m_pt->m_pa),
            convert(m_pt->m_pb),
            convert(m_pt->m_pc));
    }

    // returns a reference point which must be on the primitive
    Point reference_point() const
    { return convert(m_pt->m_pa); }
};



typedef CGAL::AABB_traits<K, My_triangle_primitive> My_AABB_traits;
typedef CGAL::AABB_tree<My_AABB_traits> Tree;

int main()
{
    My_point a(1.0, 0.0, 0.0);
    My_point b(0.0, 1.0, 0.0);
    My_point c(0.0, 0.0, 1.0);
    My_point d(0.0, 0.0, 0.0);

    std::vector<My_triangle> triangles;
    triangles.push_back(My_triangle(&a,&b,&c));
    triangles.push_back(My_triangle(&a,&b,&d));
    triangles.push_back(My_triangle(&a,&d,&c));

    // constructs AABB tree
    Tree tree(triangles.begin(),triangles.end());

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
