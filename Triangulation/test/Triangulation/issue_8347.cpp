#include <iostream>

#include <CGAL/Epick_d.h>
#include <CGAL/Epeck_d.h>
#include <CGAL/Delaunay_triangulation.h>

typedef CGAL::Dynamic_dimension_tag Dim_tag;
//typedef CGAL::Dimension_tag<2> Dim_tag;

typedef CGAL::Epick_d< Dim_tag >  Kernel;
typedef CGAL::Delaunay_triangulation<Kernel> Triangulation;

typedef typename Triangulation::Vertex_handle Vertex_handle;
typedef typename Triangulation::Point Point;

Point create_point(double x,double y) {
    std::vector<double> c;
    c.push_back(x);
    c.push_back(y);
    return Point(2,c.begin(),c.end());
}

int main() {

    Triangulation T(2);

    Point p1=create_point(2, 3);
    Point p2=create_point(6, 3);
    Point p3=create_point(0, 4);

    Vertex_handle vh1=T.insert(p1);
    T.insert(p2);
    T.remove(vh1);

    T.insert(p3);

    std::cout << "Exit normally" << std::endl;

    return 0;
}
