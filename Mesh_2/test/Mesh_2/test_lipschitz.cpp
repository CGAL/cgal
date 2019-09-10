#include "test_dependencies.h"

// and additionally
#define CGAL_SPATIAL_SEARCHING_COMMERCIAL_LICENSE 22222222
#define CGAL_APOLLONIUS_GRAPH_2_COMMERCIAL_LICENSE 22222222

#include <CGAL/config.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>

#include <CGAL/Mesh_2/Lipschitz_sizing_field_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_mesh_vertex_base_2<K>  Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K>    Fb;

typedef CGAL::Triangulation_data_structure_2<Vb, Fb>  TDS;
typedef CGAL::Exact_predicates_tag                    Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> CDT;


int main()
{
  CDT cdt;

  CGAL::Lipschitz_sizing_field_2<CDT> lip_size(cdt);
  lip_size.set_K(2.);
  CGAL_assertion(lip_size.get_K() == 2.);

}
