#include <CGAL/basic.h>
#include <iostream>
#include <CGAL/Gmpz.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Triangulation_euclidean_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>

using namespace CGAL;

typedef Homogeneous<CGAL::Gmpz>  Rp;
typedef Triangulation_euclidean_traits_xy_3<Rp>  Gt;
typedef Triangulation_vertex_base_2<Gt> Vb;
typedef Triangulation_face_base_2<Gt> Fb;
typedef Triangulation_default_data_structure_2<Gt,Vb,Fb > Tds;
typedef Delaunay_triangulation_2<Gt, Tds> Delaunay;

int main()
{
    Delaunay dt;

    Point_3<Rp> p;
    while(std::cin >> p){
      dt.insert(p);
    }
    return 0;
}

