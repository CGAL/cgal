#include "MyMink.h"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <iostream>
#include <vector>
using namespace std;

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polygon_2<Kernel> Polygon_2;
typedef CGAL::Point_2<Kernel> Point_2;
typedef CGAL::Polygon_with_holes_2<Kernel> Polygon_with_holes_2;

int main(int argc,char* argv[]){
    vector<Point_2> vertices_a;
    vertices_a.push_back(Point_2(0,0));
    vertices_a.push_back(Point_2(1,0));
    vertices_a.push_back(Point_2(1,1));
    vertices_a.push_back(Point_2(0,1));
    Polygon_2 polygon_a(vertices_a.begin(), vertices_a.end());

    vector<Point_2> vertices_b;
    vertices_b.push_back(Point_2(0,0));
    vertices_b.push_back(Point_2(1,0.5));
    vertices_b.push_back(Point_2(0,1));
    Polygon_2 polygon_b(vertices_b.begin(), vertices_b.end());

    Polygon_with_holes_2 sum = minkowski_sum_2_(polygon_a,polygon_b);
    cout << sum << endl;
}
