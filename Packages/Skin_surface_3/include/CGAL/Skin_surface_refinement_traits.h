#ifndef SKIN_SURFACE_REFINEMENT_TRAITS_H
#define SKIN_SURFACE_REFINEMENT_TRAITS_H

CGAL_BEGIN_NAMESPACE

template <class Triangulation_3, class PolyhedronTraits_3,
          class T2P_converter, class P2T_converter>
class Skin_surface_refinement_traits_3 {
public:
  typedef Triangulation_3                         Triangulation;
  typedef PolyhedronTraits_3                      Polyhedron_traits;

  typedef typename Polyhedron_traits::Point_3     Polyhedron_point;
  typedef typename Polyhedron_traits::Segment_3   Polyhedron_segment;
  typedef typename Polyhedron_traits::Line_3      Polyhedron_line;
  typedef typename Polyhedron_traits::Plane_3     Polyhedron_plane;

  typedef typename Triangulation::Cell_handle     Triang_cell_handle;

  Skin_surface_refinement_traits_3(
    Triangulation const& t, 
    T2P_converter t2p_converter, 
    P2T_converter p2t_converter) 
  : t(t), t2p_converter(t2p_converter), p2t_converter(p2t_converter) {
    
  }
    
  // Additional functions, not belonging to the traits concept:
  //Polyhedron_rt value(const Cell_handle &ch, const Polyhedron_point &p) const {
    //return ch->surf->value(p);
  //}
//  Polyhedron_rt value(const Cell_handle &ch, const Triang_Point &p) const {
//    return ch->surf->value(converter(p));
  //}
  Polyhedron_point to_surface_along_transversal_segment(
    Polyhedron_point const &p, Triang_cell_handle ch) {
      
    Polyhedron_segment transversal_segment = get_transversal_segment(ch,p);
      
    return ch->surf->to_surface(
      transversal_segment.begin(), transversal_segment.end());
  }
  std::pair<Polyhedron_point, Polyhedron_point> get_transversal_segment(
    Triang_cell_handle ch, Polyhedron_point const &p) {
    // Compute signs on vertices and sort them:
    int nIn = 0;
    for (int i=0; i<4; i++) {
      if (is_inside(cit,i)) {
        sortedV[nIn] = i; nIn++;
      } else {
        sortedV[3-i+nIn] = i;
      }
    }
    Polyhedron_point begin_point, end_point;
    if (nIn==1) {
      begin_point = P2T_converter(ch->vertex(sortedV[0])),
      end_point = intersection(
        Polyhedron_plane(
          P2T_converter(ch->vertex(sortedV[1])),
          P2T_converter(ch->vertex(sortedV[2])),
          P2T_converter(ch->vertex(sortedV[3]))),
        Polyhedron_line(
          begin_point,
          P2T_converter(ch->vertex(p))));
    } else if (nIn==2) {
      begin_point = intersection(
        Polyhedron_plane(
          P2T_converter(ch->vertex(sortedV[2])),
          P2T_converter(ch->vertex(sortedV[3])),
          p),
        Polyhedron_line(
          P2T_converter(ch->vertex(sortedV[0])),
          P2T_converter(ch->vertex(sortedV[1]))));
      begin_point = intersection(
        Polyhedron_plane(
          P2T_converter(ch->vertex(sortedV[0])),
          P2T_converter(ch->vertex(sortedV[1])),
          p),
        Polyhedron_line(
          P2T_converter(ch->vertex(sortedV[2])),
          P2T_converter(ch->vertex(sortedV[3]))));
    } else if (nIn==3) {
      end_point = P2T_converter(ch->vertex(sortedV[3])),
      begin_point = intersection(
        Polyhedron_plane(
          P2T_converter(ch->vertex(sortedV[0])),
          P2T_converter(ch->vertex(sortedV[1])),
          P2T_converter(ch->vertex(sortedV[2]))),
        Polyhedron_line(
          end_point,
          P2T_converter(ch->vertex(p))));
    } else {
      CGAL_assertion(false);
    }
    return Polyhedron_segment(begin_point, end_point);
  }
  
  Polyhedron_point 
  intersection(Polyhedron_plane const &plane, Polyhedron_line const &line) {
    Polyhedron_point p;

    Object result = CGALi::intersection(plane, line);
    if ( !CGAL::assign(p, result) )
      CGAL_assertion_msg(false,"intersection: no intersection.");
    return p;
  }
  Sign sign(Triang_cell_handle ch, int i) const {
    return CGAL_NTS sign(
       ch->surf->value(converter(ch->vertex(i)->point())) - iso_value);
  }  
private:
  Triangulation const &t;
  T2P_converter const &t2p_converter;
  P2T_converter const &p2t_converter;
};

CGAL_END_NAMESPACE

#endif // SKIN_SURFACE_REFINEMENT_TRAITS_H
