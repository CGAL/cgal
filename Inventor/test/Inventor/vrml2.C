// file: test/Inventor/VRML2.C

#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/VRML_2_ostream.h>
#include <iostream>

typedef CGAL::Simple_cartesian<double>  Kernel;
typedef CGAL::Point_3<Kernel>           Point;
typedef CGAL::Segment_3<Kernel>         Segment;
typedef CGAL::Triangle_3<Kernel>        Triangle;
typedef CGAL::Tetrahedron_3<Kernel>     Tetrahedron;
typedef CGAL::Sphere_3<Kernel>          Sphere;

int main() {
    Point p1(10, 15, 0), p2(0, 0, 0), 
      p3(-10, 0, -15), p4(15, 15, 15);
    Segment s1(p1, p2);
    Triangle t1(p1, p2, p3);
    Tetrahedron tet1(p1, p2, p3, p4);
    Sphere s(p1, p2);
    CGAL::VRML_2_ostream out( std::cout);
    out << s1 << t1 << p1 << tet1 << s;
    return 0;
}
