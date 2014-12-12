#include <CGAL/Exact_integer.h>
#include<CGAL/Homogeneous.h>
#include<CGAL/Nef_polyhedron_3.h>
#include<CGAL/IO/Nef_polyhedron_iostream_3.h>
#include<CGAL/Polyhedron_3.h>
#include<CGAL/convex_decomposition_3.h>
#include<CGAL/Nef_nary_union_3.h>
#include<CGAL/Nef_3/SNC_indexed_items.h>
#include<CGAL/Nef_S2/SM_const_decorator.h>
#include<CGAL/Minkowski_sum_3/Gaussian_map.h>
#include<CGAL/Minkowski_sum_3/Gaussian_map_to_nef_3.h>
#include<fstream>

typedef CGAL::Exact_integer RT;
typedef CGAL::Homogeneous<RT> Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel, CGAL::SNC_indexed_items> Nef_polyhedron_3;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef CGAL::Nef_nary_union_3<Nef_polyhedron_3> Nary_union;
typedef Nef_polyhedron_3::SFace_const_iterator SFace_const_iterator;
typedef Nef_polyhedron_3::SNC_structure SNC_structure;
typedef SNC_structure::Sphere_map Sphere_map;
typedef CGAL::SM_const_decorator<Sphere_map> SM_const_decorator;

void test_convex_parts(Nef_polyhedron_3& N)
{
  //  Nef_polyhedron_3 N(C); //TODO: test
  CGAL::convex_decomposition_3(N);
  std::list<Nef_polyhedron_3> convex_parts;
  Nef_polyhedron_3::Volume_const_iterator ci;
  for(ci = ++N.volumes_begin(); ci != N.volumes_end(); ++ci) {
    if(ci->mark()) {
      Polyhedron_3 P;
      N.convert_inner_shell_to_polyhedron(ci->shells_begin(), P);
      Nef_polyhedron_3 tmp0(P), tmp1;

      CGAL::Gaussian_map<Kernel,Nef_polyhedron_3> G(N, ci);
      CGAL::Gaussian_map_to_nef_3<Nef_polyhedron_3> 
	Convertor(G);
      tmp1.delegate(Convertor, true);

      assert(tmp1.is_valid());      
      assert(tmp1.closure().symmetric_difference(tmp0).is_empty());
      delete G.sphere_map();
    }
  }

  Nef_polyhedron_3::Vertex_const_iterator vi;
  for(vi = N.vertices_begin(); vi != N.vertices_end(); ++vi) {
    if(vi->number_of_sfaces() == 1 &&
       vi->number_of_svertices() == 0 &&
       vi->mark() &&
       !vi->sfaces_begin()->mark()) {
      Nef_polyhedron_3 tmp;
      CGAL::Gaussian_map<Kernel,Nef_polyhedron_3> G(vi);
      CGAL::Gaussian_map_to_nef_3<Nef_polyhedron_3> 
	Convertor(G);
      tmp.delegate(Convertor, true);
      assert(tmp.is_valid());
      assert(tmp.number_of_vertices() == 1);
      assert(tmp.number_of_volumes() == 1);
      assert(tmp.vertices_begin()->mark());
      assert(!tmp.volumes_begin()->mark());
      assert(tmp.vertices_begin()->point() == vi->point());

      delete G.sphere_map();

      std::cerr << "single vertex " << std::endl;
    }
  }

  Nef_polyhedron_3::Halfedge_const_iterator eci;
  for(eci = N.halfedges_begin(); eci != N.halfedges_end(); ++eci) {
    if(eci->is_twin()) continue;
    SM_const_decorator SD(&*eci->source());
    if(!SD.is_isolated(eci)) continue;
    Nef_polyhedron_3 tmp;
    CGAL::Gaussian_map<Kernel,Nef_polyhedron_3> G(eci);
    CGAL::Gaussian_map_to_nef_3<Nef_polyhedron_3> 
      Convertor(G);
    tmp.delegate(Convertor, true);
    assert(tmp.is_valid());
    assert(tmp.number_of_vertices() == 2);
    assert(tmp.number_of_halfedges() == 2);
    assert(tmp.number_of_volumes() == 1);
    assert(tmp.vertices_begin()->mark());
    assert((++tmp.vertices_begin())->mark());
    assert(!tmp.volumes_begin()->mark());
    assert(tmp.vertices_begin()->point() ==
		   eci->source()->point() ||
		   tmp.vertices_begin()->point() ==
		   eci->twin()->source()->point());
    delete G.sphere_map();

    std::cerr << "single edge " << std::endl;
  }
  
  Nef_polyhedron_3::Halffacet_const_iterator fci;
  for(fci = N.halffacets_begin(); fci != N.halffacets_end(); ++fci) {
    if(fci->is_twin()) continue;  
    if(fci->incident_volume() !=
       fci->twin()->incident_volume()) continue;
    Nef_polyhedron_3 tmp;
    CGAL::Gaussian_map<Kernel, Nef_polyhedron_3> G(fci);
    CGAL::Gaussian_map_to_nef_3<Nef_polyhedron_3> 
      Convertor(G);
    tmp.delegate(Convertor, true);
    assert(tmp.is_valid());
    assert(tmp.number_of_halffacets() == 2);
    assert(tmp.number_of_volumes() == 1);
    assert(tmp.number_of_vertices()*2 == tmp.number_of_halfedges());
    assert(tmp.halffacets_begin()->mark());
    assert(!tmp.volumes_begin()->mark());
    assert(tmp.halffacets_begin()->plane() == fci->plane() ||
		   tmp.halffacets_begin()->plane() == fci->twin()->plane());
    delete G.sphere_map();

    std::cerr << "single facet " << std::endl;
  }
}

int main()
{
  Nef_polyhedron_3 N;

  std::ifstream in0("star.nef3");
  in0 >> N;
  test_convex_parts(N);
  
  std::ifstream in1("single_vertex.nef3");
  in1 >> N;
  test_convex_parts(N);

  std::ifstream in2("single_edge.nef3");
  in2 >> N;
  test_convex_parts(N);

  std::ifstream in3("single_facet10.nef3");
  in3 >> N;
  test_convex_parts(N);
}
