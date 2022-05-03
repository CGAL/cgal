#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/boost/graph/generators.h>


#include <iostream>
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<K::Point_3> Polygon_mesh;

typedef CGAL::Nef_polyhedron_3<K> Nef_polyhedron;

int main(int /* argc */, char** /* argv[] */)
{
  Point_3 p0(1,1,1), p1(2,1,1), p2(2,2,1), p3(1,2,1), p4(1,2,2), p5(1,1,2), p6(2,1,2), p7(2,2,2);
  Point_3 q0(0,0,0), q1(3,0,0), q2(3,3,0), q3(0,3,0), q4(0,3,3), q5(0,0,3), q6(3,0,3), q7(3,3,3);

    Polygon_mesh A1, A2;

    make_hexahedron(p0, p1, p2, p3, p4, p5, p6, p7, A1);
    make_hexahedron(q0, q1, q2, q3, q4, q5, q6, q7, A2);
    Nef_polyhedron a1(A1), a2(A2);

    assert(a1 < a2);
    assert(a2 > a1);

    return 0;
}
