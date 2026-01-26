#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Conforming_constrained_Delaunay_triangulation_3.h>
#include <CGAL/property_map.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<std::size_t, K> Vb;
typedef CGAL::Conforming_constrained_Delaunay_triangulation_vertex_base_3<K, Vb> CVb;
typedef CGAL::Triangulation_cell_base_with_info_3<bool, K> Cb;
typedef CGAL::Conforming_constrained_Delaunay_triangulation_cell_base_3<K, Cb>  CCb;

typedef CGAL::Triangulation_data_structure_3<CVb, CCb> Tds;
typedef CGAL::Triangulation_3<K, Tds> Trig_3;
typedef CGAL::Conforming_constrained_Delaunay_triangulation_3<K, Trig_3> CstrDelaunayTrig_3;

int main()
{
    std::vector<std::pair<K::Point_3, std::size_t>> points;

    std::vector<std::vector<std::size_t>> polygons;

    using First_pmap = CGAL::First_of_pair_property_map<std::pair<K::Point_3, std::size_t>>;
    //Initialize the point list
    points.push_back({K::Point_3(0, 0, 0), 1});
    points.push_back({K::Point_3(1, 0, 0), 2});
    points.push_back({K::Point_3(0, 1, 0), 3});
    points.push_back({K::Point_3(0, 0, 1), 4});
    polygons.push_back({0, 1, 2});
    polygons.push_back({0, 1, 3});
    CstrDelaunayTrig_3 cdt(points, polygons, CGAL::parameters::point_map(First_pmap()));
    for(const auto& v : cdt.triangulation().finite_vertex_handles())
    {
        std::cout << "Vertex: " << v->point() << ", Info: " << v->info() << std::endl;
    }

    return 0;
}
