// Mesh_3 bug: The surface of c3t3 has holes
//
//   https://github.com/CGAL/cgal/issues/1554
//
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>

#include <fstream>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;

typedef CGAL::Parallel_if_available_tag Concurrency_tag;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,
                                   CGAL::Default,
                                   Concurrency_tag>::type  Tr;
typedef Tr::Geom_traits                                    Gt;

typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr,
                                                Mesh_domain::Corner_index,
                                                Mesh_domain::Curve_index> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

template < class ConstructWCircumcenter >
Tr::Bare_point wc_circumcenter(const Tr& tr,
                               Tr::Cell_handle ch,
                               ConstructWCircumcenter construct_w_circumcenter,
                               bool force_exact = false)
{
  return construct_w_circumcenter(tr.point(ch, 0), tr.point(ch, 1),
                                  tr.point(ch, 2), tr.point(ch, 3), force_exact);
}

int main(int argc, char*argv[])
{
  const char* fname = (argc>1)?argv[1]:"data/fandisk.off";
  // Create domain
  std::ifstream in(fname);
  Polyhedron poly;
  in >> poly;
  Mesh_domain domain(poly);

  // Get sharp features
  domain.detect_features();

  // Mesh criteria
  Mesh_criteria criteria(edge_size = 0.037,
                         facet_angle = 25,
                         facet_size = 0.037,
                         facet_distance = 0.0037,
                         cell_radius_edge_ratio = 3);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());

  Gt::Construct_weighted_circumcenter_3 w_circumcenter =
      c3t3.triangulation().geom_traits().construct_weighted_circumcenter_3_object();

  int return_code = 0;
  for(C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin();
      cit != c3t3.cells_in_complex_end();
      ++cit)
  {
    for(int i = 0; i < 4; ++i) {
      typedef Tr::Facet Facet;
      typedef Tr::Cell_handle Cell_handle;
      Cell_handle ch = cit;
      Cell_handle nch = ch->neighbor(i);
      if(c3t3.is_in_complex(nch)) continue;
      if(!c3t3.is_in_complex(Facet(ch, i))) {
        return_code = 1;
        std::cout << "ERROR with facet: 3 \n";
        K::Weighted_point_3 p[3];
        for(int j = 0; j < 3; ++j) {
          p[j] = c3t3.triangulation().point(ch, c3t3.triangulation().vertex_triple_index(i, j));
          std::cout << "  " << p[j] << "\n";
        }
        int n_i = nch->index(ch);
        K::Weighted_point_3 a = c3t3.triangulation().point(ch, i);
        K::Weighted_point_3 b = c3t3.triangulation().point(nch, n_i);
        std::cout << a << "\n";
        std::cout << b << "\n";
        std::cout << "    " << wc_circumcenter(c3t3.triangulation(), ch, w_circumcenter, true) << std::endl;
        std::cout << "    " << wc_circumcenter(c3t3.triangulation(), nch, w_circumcenter, true) << std::endl;
      }
    }
  }
  std::cout << c3t3.triangulation().number_of_vertices() << std::endl;
  // Output
  return return_code;
}
