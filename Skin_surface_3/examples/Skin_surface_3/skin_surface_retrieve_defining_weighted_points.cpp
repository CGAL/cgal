#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Skin_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_skin_surface_3.h>
#include <list>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Skin_surface_traits_3<K> Traits;
typedef CGAL::Skin_surface_3<Traits> Skin_surface_3;
typedef Skin_surface_3::FT FT;
typedef Skin_surface_3::Weighted_point Weighted_point;
typedef Weighted_point::Point Bare_point;
typedef CGAL::Skin_surface_polyhedral_items_3<Skin_surface_3> Polyhedral_items;
typedef CGAL::Polyhedron_3<K, Polyhedral_items> Polyhedron;

int main() {
  std::list<Weighted_point> l;
  FT shrinkfactor = 0.5;

  l.push_front(Weighted_point(Bare_point(0, 0, 0), 1));
  l.push_front(Weighted_point(Bare_point(0, 1, 0), 2));
  l.push_front(Weighted_point(Bare_point(0, 0, 2), 1));

  Skin_surface_3 skin_surface(l.begin(), l.end(), shrinkfactor);

  Polyhedron p;
  CGAL::mesh_skin_surface_3(skin_surface, p);

  for (Polyhedron::Facet_iterator fit = p.facets_begin(); fit != p.facets_end(); ++fit) {

    // Retrieve the generating simplex from the regular triangulation
    const Skin_surface_3::Simplex &s = fit->containing_simplex();

    // Get the weighted points from the simplex
    Skin_surface_3::Weighted_point wp[4];
    switch (s.dimension()) {
    case 0: {
      Skin_surface_3::Vertex_handle vh = s;
      wp[0] = vh->point();
      break;
    }
    case 1: {
      Skin_surface_3::Edge e = s;
      wp[0] = e.first->vertex(e.second)->point();
      wp[1] = e.first->vertex(e.third)->point();
      break;
    }
    case 2: {
      Skin_surface_3::Facet f = s;
      wp[0] = f.first->vertex((f.second + 1) & 3)->point();
      wp[1] = f.first->vertex((f.second + 2) & 3)->point();
      wp[2] = f.first->vertex((f.second + 3) & 3)->point();
      break;
    }
    case 3: {
      Skin_surface_3::Cell_handle ch = s;
      wp[0] = ch->vertex(0)->point();
      wp[1] = ch->vertex(1)->point();
      wp[2] = ch->vertex(2)->point();
      wp[3] = ch->vertex(3)->point();
      break;
    }
    }
  }

  return 0;
}
