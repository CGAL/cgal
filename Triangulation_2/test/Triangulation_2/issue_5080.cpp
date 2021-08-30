// clang++ a01.cpp -o a01 -I/home/jherring/code/dev/third-party/current/include -lgmp -lmpfr

#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

template<typename VertexInfoType, typename FaceInfoType, typename Itag>
struct TriangulationTraitsTemplate_2 {

    using K = CGAL::Exact_predicates_inexact_constructions_kernel;

    using Vbb = CGAL::Triangulation_vertex_base_with_info_2<VertexInfoType, K>;
    using Vb = CGAL::Triangulation_hierarchy_vertex_base_2<Vbb>;

    using Fbb = CGAL::Triangulation_face_base_with_info_2<FaceInfoType, K>;
    using Fb = CGAL::Constrained_triangulation_face_base_2<K, Fbb>;

    using TDS = CGAL::Triangulation_data_structure_2<Vb, Fb>;

    using CDTBase = CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>;

    using CDTHierarchy = CGAL::Triangulation_hierarchy_2<CDTBase>;

    using CDT = CGAL::Constrained_triangulation_plus_2<CDTHierarchy>;

};

template<typename VertexInfoType, typename FaceInfoType>
struct TriangulationTraitsTemplate_2_AllowIntersections {
    using CDT = typename TriangulationTraitsTemplate_2<
        VertexInfoType,
        FaceInfoType,
        CGAL::Exact_predicates_tag>::CDT;
};

template<int DIM, typename VertexInfoType, typename FaceInfoType>
struct TriangulationTraitsTemplate {};

template<typename VertexInfoType, typename FaceInfoType>
struct TriangulationTraitsTemplate<2, VertexInfoType, FaceInfoType>{
    using CDT = typename TriangulationTraitsTemplate_2_AllowIntersections<
        VertexInfoType,
        FaceInfoType>::CDT;
};

struct VertexInfo
{};

struct FaceInfo
{};

using Triangulation = TriangulationTraitsTemplate<2, VertexInfo, FaceInfo>::CDT;
using Point = Triangulation::Point;
using Edge = Triangulation::Edge;

int main()
{
    Triangulation triangulation;

    Point p0{ -1.5617169965001162502e-21, -0.0059749999999959748503 };
    std::vector<std::array<Point, 2>> as = {
        { p0, Point{ 1.8241074394927896017e-06, -0.0059800116939744598493 } },
        { p0, Point{ 1.8241073772900120186e-06, -0.0059699883060029008269 } },
        { p0, Point{ 5.0116939794452514324e-06, -0.00597682410742579346   } },
        { p0, Point{ 5.0116939818259595034e-06, -0.0059731758925807469998 } }
    };
    std::vector<std::array<Point, 2>> bs = {
        { Point{  0.0, -0.0060000000000000001249 }, Point{  0.0, -0.00596000000000000002 } }
    };

    for (const auto& a: as) {
        auto v1 = triangulation.insert(a[0]);
        auto v2 = triangulation.insert(a[1]);
        triangulation.insert_constraint(v1, v2);
    }

    for (const auto& b: bs) {
        auto v1 = triangulation.insert(b[0]);
        auto v2 = triangulation.insert(b[1]);
        triangulation.insert_constraint(v1, v2);
    }

    return 0;
}
