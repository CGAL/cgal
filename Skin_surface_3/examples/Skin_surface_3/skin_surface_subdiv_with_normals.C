// examples/Skin_surface_3/skin_surface_subdiv_with_normals.C
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Skin_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Skin_surface_polyhedral_items_3.h>
#include <CGAL/mesh_skin_surface_3.h>
#include <CGAL/subdivide_skin_surface_mesh_3.h>
#include <list>

#include <fstream>
#include <CGAL/IO/Polyhedron_iostream.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_euclidean_traits_3<K>   Traits;
typedef CGAL::Skin_surface_3<Traits>                        Skin_surface_3;
typedef Skin_surface_3::FT                                  FT;
typedef Skin_surface_3::Weighted_point                      Weighted_point;
typedef Skin_surface_3::Bare_point                          Bare_point;

// Each facet has a pointer to the tetrahedron of the TMC it is contained in
typedef Skin_surface_3::Triangulated_mixed_complex          TMC;
typedef CGAL::Skin_surface_polyhedral_items_3<TMC>          Poly_items;
typedef CGAL::Polyhedron_3<K,Poly_items>                    Polyhedron;


template <class Polyhedron, class SkinSurface>
/// Write polyhedron with normals:
void write_polyhedron_with_normals(Polyhedron &p, 
				   SkinSurface &skin, 
				   std::ostream &out)
{
  typedef typename Polyhedron::Vertex_iterator                  Vertex_iterator;
  typedef typename Polyhedron::Facet_iterator                   Facet_iterator;
  typedef typename Polyhedron::Halfedge_around_facet_circulator HFC;
  typedef typename Polyhedron::Vertex_handle                    Vertex_handle;

  // Write header
  out << "NOFF " << p.size_of_vertices ()
      << " " << p.size_of_facets()
      << " " << p.size_of_halfedges()
      << std::endl;

  // Write vertices
  typedef CGAL::Skin_surface_subdivision_policy_base_3<Polyhedron, SkinSurface> 
    Subdivision_policy;
  Subdivision_policy *policy = get_subdivision_policy(p, skin);
  for (Vertex_iterator vit = p.vertices_begin();
       vit != p.vertices_end(); vit ++) {
    out << vit->point() << " "
 	<< policy->normal(vit)
	<< std::endl;
  }

  // Write faces
  CGAL::Inverse_index<Vertex_handle> index(p.vertices_begin(),
					   p.vertices_end());
  for(Facet_iterator fi = p.facets_begin();
      fi != p.facets_end(); ++fi) {
    HFC hc = fi->facet_begin();
    HFC hc_end = hc;
    std::size_t n = circulator_size( hc);
    out << n;
    do {
      Vertex_handle vh = (*hc).vertex();
      out << " " << index[vh];
    } while (++hc != hc_end);
    out << "\n";
  }
}  

int main(int argc, char *argv[]) {
  std::list<Weighted_point> l;
  FT                        shrinkfactor = 0.5;
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <infile>" << std::endl;
    exit(0);
  }

  std::ifstream in(argv[1]);
  Weighted_point wp;
  while (in >> wp) l.push_front(wp);

  Skin_surface_3 skin_surface(l.begin(), l.end(), shrinkfactor);

  Polyhedron p;
  CGAL::mesh_skin_surface_3(skin_surface, p);

  CGAL::subdivide_skin_surface_mesh_3(p, skin_surface);
  
  std::ofstream out("mesh.off");
  write_polyhedron_with_normals(p, skin_surface, out);

  return 0;
}
