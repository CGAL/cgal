// file : example/Triangulation_2/terrain.C

#include <CGAL/Gmpz.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Triangulation_euclidean_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>

typedef CGAL::Homogeneous<CGAL::Gmpz>  Rp;
typedef CGAL::Triangulation_euclidean_traits_xy_3<Rp>  Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Triangulation_face_base_2<Gt> Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb > Tds;
typedef CGAL::Delaunay_triangulation_2<Gt, Tds> Delaunay;

int main()
{
    Delaunay dt;

    Rp::Point_3  p;
    while(std::cin >> p){
      dt.insert(p);
    }
    return 0;
}

