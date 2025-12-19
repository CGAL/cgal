template <typename Vertex_handle>
auto display_vert(Vertex_handle v) {
  std::stringstream os;
  os.precision(17);
  if(v->time_stamp() == 0) {
    os << "inf";
  } else {
    os << '#' << v->time_stamp() << "=(" << v->point() << ")";
  }
  return os.str();
};

template <typename DT>
struct Debug_simplex {
  using Cell_handle = typename DT::Cell_handle;
  using Edge = typename DT::Edge;
  using Facet = typename DT::Facet;
  using Vertex_handle = typename DT::Vertex_handle;
  using Simplex = typename DT::Simplex;

  Simplex simplex;

  template<typename Dt, typename CharT, typename Traits>
  friend
  std::basic_ostream<CharT, Traits>&
  operator<<(std::basic_ostream<CharT, Traits>& os, const Debug_simplex<Dt>& d) {
    auto&& simplex = d.simplex;
    switch(simplex.dimension()) {
      case 0: {
        os << "- vertex " << display_vert(static_cast<Vertex_handle>(simplex));
        break;
      }
      case 1: {
        const auto [c, index1, index2] = static_cast<Edge>(simplex);
        os << "- edge "
           << display_vert(c->vertex(index1)) << " - "
           << display_vert(c->vertex(index2));
        break;
      }
      case 2: {
        const auto [c, index] = static_cast<Facet>(simplex);
        os << "- facet "
           << display_vert(c->vertex(DT::vertex_triple_index(index, 0))) << " - "
           << display_vert(c->vertex(DT::vertex_triple_index(index, 1))) << " - "
           << display_vert(c->vertex(DT::vertex_triple_index(index, 2)));
        break;
      }
      case 3: {
        const auto c = static_cast<Cell_handle>(simplex);
        os << "- cell "
           << display_vert(c->vertex(0)) << " - "
           << display_vert(c->vertex(1)) << " - "
           << display_vert(c->vertex(2)) << " - "
           << display_vert(c->vertex(3));
        break;
      }
      default: CGAL_assume(false);
    }
    return os;
  };
};

#include <CGAL/Triangulation_simplex_3.h>

template <typename Triangulation>
auto debug_simplex(CGAL::Triangulation_simplex_3<Triangulation> simplex) {
  return Debug_simplex<Triangulation>{simplex};
}
