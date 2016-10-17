#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>

// Domain 
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;

#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;

typedef CGAL::Mesh_complex_3_in_triangulation_3<
  Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_segment_index> C3t3;
typedef Tr::Geom_traits Gt;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

Gt::Construct_weighted_circumcenter_3 w_circumcenter = Gt().construct_weighted_circumcenter_3_object();

K::Point_3 wc_circumcenter(Tr::Cell_handle ch, bool force_exact = false) {
  return w_circumcenter(ch->vertex(0)->point(),
                        ch->vertex(1)->point(),
                        ch->vertex(2)->point(),
                        ch->vertex(3)->point(), force_exact
                        );
}

int main(int argc, char*argv[])
{
  const char* fname = (argc>1)?argv[1]:"data/fandisk.off";
  // Create domain
  Mesh_domain domain(fname);
  
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
        K::Point_3 p[3];
        for(int j = 0; j < 3; ++j) {
          p[j] = ch->vertex(c3t3.triangulation().vertex_triple_index(i, j))->point().point();
          std::cout << "  " << p[j] << "\n";
        }
        int n_i = nch->index(ch);
        K::Point_3 a = ch->vertex(i)->point().point();
        K::Point_3 b = nch->vertex(n_i)->point().point();
        std::cout << a << "\n";
        std::cout << b << "\n";
        std::cout << "    " << wc_circumcenter(ch, true) << std::endl;
        std::cout << "    " << wc_circumcenter(nch, true) << std::endl;
      }
    }
  }
  std::cout << c3t3.triangulation().number_of_vertices() << std::endl;
  // Output
  return return_code;
}
