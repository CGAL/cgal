#ifndef SKIN_SURFACE_REFINEMENT_TRAITS_H
#define SKIN_SURFACE_REFINEMENT_TRAITS_H

template <class Triangulation_3, class PolyhedronTraits_3,
          class T2P_converter, class P2T_converter>
Skin_surface_refinement_traits_3 {
public:
  typedef Triangulation_3                         Triangulation;
  typedef PolyhedronTraits_3                      Polyhedron_traits;

  Skin_surface_refinement_traits_3(Triangulation const& t) : t(t) {
    
  }
    
  // Additional functions, not belonging to the traits concept:
  HDS_rt value(const Cell_handle &ch, const HDS_point &p) const {
    return ch->surf->value(p);
  }
  HDS_rt value(const Cell_handle &ch, const Triang_Point &p) const {
    return ch->surf->value(converter(p));
  }
  HDS_point to_surface_along_transversal_segment(
    HDS_point const &p, Cell_handle ch) {
    std::pair<HDS_point, HDS_point>
  }
  std::pair<HDS_point, HDS_point> get_transversal_segment(
    Cell_handle ch, HDS_point const &p) {
  }
private:
  Triangulation const &t;
};

#endif // SKIN_SURFACE_REFINEMENT_TRAITS_H
