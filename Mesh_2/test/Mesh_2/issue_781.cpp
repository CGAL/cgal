#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

// mesh refinement
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel                                 meshTriKernel;
typedef meshTriKernel::Point_2                                                              meshTriPoint;

typedef CGAL::Triangulation_vertex_base_2<meshTriKernel>                                    meshTriVertexBase;

typedef CGAL::Constrained_triangulation_face_base_2<meshTriKernel>                          Fbb;
typedef CGAL::Delaunay_mesh_face_base_2<meshTriKernel,Fbb>                                  meshTriFaceBase;

typedef CGAL::Triangulation_data_structure_2<meshTriVertexBase,meshTriFaceBase>             meshTriTDS;
typedef CGAL::Exact_intersections_tag                                                       meshTriItag;
typedef CGAL::Constrained_Delaunay_triangulation_2<meshTriKernel, meshTriTDS, meshTriItag>  meshTriCDT;

typedef CGAL::Delaunay_mesh_size_criteria_2<meshTriCDT>                                 meshCriteria;
typedef CGAL::Delaunay_mesher_2<meshTriCDT, meshCriteria>                               meshRefiner;

int main(int argc, char* argv[])
{
    std::cerr.precision(17);
    meshTriCDT cdt;

    meshTriPoint pt(5.4691594172333904, 44.256641611715409);

    cdt.insert( pt );

    meshTriPoint t(5.4693788499999929, 44.256578099999999);

    meshTriPoint m(5.4691178249999917, 44.256653649999997);

    cdt.insert_constraint( m,t );

    assert(cdt.is_valid());
    meshRefiner mesher(cdt);
    mesher.set_criteria(meshCriteria(0.125));

    std::cout << "refine mesh" << std::endl;
    mesher.refine_mesh();
    std::cout << "complete" << std::endl;

    return 0;
}
