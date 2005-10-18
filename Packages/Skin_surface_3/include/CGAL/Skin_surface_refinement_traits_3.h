#ifndef SKIN_SURFACE_REFINEMENT_TRAITS_H
#define SKIN_SURFACE_REFINEMENT_TRAITS_H

#include <CGAL/intersection_3_1.h>

CGAL_BEGIN_NAMESPACE

template <class Triangulation_3, class PolyhedronTraits_3,
          class T2P_converter, class P2T_converter>
class Skin_surface_refinement_traits_3 {
public:
  typedef Triangulation_3                         Triangulation;
  typedef PolyhedronTraits_3                      Polyhedron_traits;

  typedef typename Polyhedron_traits::RT          Polyhedron_rt;
  typedef typename Polyhedron_traits::Point_3     Polyhedron_point;
  typedef typename Polyhedron_traits::Segment_3   Polyhedron_segment;
  typedef typename Polyhedron_traits::Line_3      Polyhedron_line;
  typedef typename Polyhedron_traits::Vector_3    Polyhedron_vector;
  typedef typename Polyhedron_traits::Plane_3     Polyhedron_plane;

  typedef typename Triangulation::Cell_handle     Triang_cell_handle;
  typedef typename Triangulation::Vertex_handle   Triang_vertex_handle;

  typedef std::map<Triang_vertex_handle,bool>     Triang_vertex_map;
  typedef typename Triang_vertex_map::iterator    Triang_vertex_map_it;


  Skin_surface_refinement_traits_3(
    Triangulation const& t, 
    Polyhedron_rt iso=0,
    T2P_converter t2p_converter = T2P_converter(), 
    P2T_converter p2t_converter = P2T_converter())
  : t(t), iso_value(iso), t2p_converter(t2p_converter), p2t_converter(p2t_converter) {
    
  }
    
  Polyhedron_point to_surface_along_transversal_segment(
    Polyhedron_point &p) {
    Triang_cell_handle ch = t.locate(p2t_converter(p));
    return to_surface_along_transversal_segment(p,ch);
  }
  Polyhedron_point to_surface_along_transversal_segment(
    Polyhedron_point &p, Triang_cell_handle ch) {
      
    Polyhedron_segment transversal_segment = get_transversal_segment(ch,p);
      
    return ch->surf->to_surface(
      transversal_segment.start(), transversal_segment.end());
  }
  Polyhedron_vector normal(Polyhedron_point const &p) {
    Triang_cell_handle ch = t.locate(p2t_converter(p));
    return ch->surf->normal(p);
  }
  int dimension(Polyhedron_point const &p) {
    Triang_cell_handle ch = t.locate(p2t_converter(p));
    return ch->surf->dimension();
  }
private:
  Polyhedron_segment get_transversal_segment(
    Triang_cell_handle ch, Polyhedron_point &p) {
    // Compute signs on vertices and sort them:
    int nIn = 0;
    int sortedV[4];
    for (int i=0; i<4; i++) {
      if (is_inside(ch,i)) {
        sortedV[nIn] = i; nIn++;
      } else {
        sortedV[3-i+nIn] = i;
      }
    }
    Polyhedron_point begin_point, end_point;
    Object obj;
    if (nIn==1) {
      begin_point = t2p_converter(ch->vertex(sortedV[0])->point());
      obj = CGAL::intersection(
        Polyhedron_plane(
          t2p_converter(ch->vertex(sortedV[1])->point()),
          t2p_converter(ch->vertex(sortedV[2])->point()),
          t2p_converter(ch->vertex(sortedV[3])->point())),
        Polyhedron_line(begin_point, p));
      if ( !assign(end_point, obj) )
        CGAL_assertion_msg(false,"intersection: no intersection.");
    } else if (nIn==2) {
      obj = CGAL::intersection(
        Polyhedron_plane(
          t2p_converter(ch->vertex(sortedV[2])->point()),
          t2p_converter(ch->vertex(sortedV[3])->point()),
          p),
        Polyhedron_line(
          t2p_converter(ch->vertex(sortedV[0])->point()),
          t2p_converter(ch->vertex(sortedV[1])->point())));
      if ( !assign(begin_point, obj) )
        CGAL_assertion_msg(false,"intersection: no intersection.");
      obj = CGAL::intersection(
        Polyhedron_plane(
          t2p_converter(ch->vertex(sortedV[0])->point()),
          t2p_converter(ch->vertex(sortedV[1])->point()),
          p),
        Polyhedron_line(
          t2p_converter(ch->vertex(sortedV[2])->point()),
          t2p_converter(ch->vertex(sortedV[3])->point())));
      if ( !assign(end_point, obj) )
        CGAL_assertion_msg(false,"intersection: no intersection.");
    } else if (nIn==3) {
      end_point = t2p_converter(ch->vertex(sortedV[3])->point());
      obj = CGAL::intersection(
        Polyhedron_plane(
          t2p_converter(ch->vertex(sortedV[0])->point()),
          t2p_converter(ch->vertex(sortedV[1])->point()),
          t2p_converter(ch->vertex(sortedV[2])->point())),
        Polyhedron_line(end_point, p));
      if ( !assign(begin_point, obj) )
        CGAL_assertion_msg(false,"intersection: no intersection.");
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

  Sign sign(Triang_cell_handle ch, int i) {
    return CGAL_NTS sign(
       ch->surf->value(t2p_converter(ch->vertex(i)->point())) - iso_value);
  }  
  bool is_inside(Triang_cell_handle ch, int i) {
    //return (sign(ch,i) == POSITIVE);
    Triang_vertex_map_it it = triang_vertex_signs.find(ch->vertex(i));
    
    if (it == triang_vertex_signs.end()) {
      bool side = (sign(ch,i) == POSITIVE);
      triang_vertex_signs[ch->vertex(i)] = side;
      CGAL_assertion(triang_vertex_signs[ch->vertex(i)] == side);
      return side;
    } else {
      return it->second;
    }
  }

  
private:
  Triangulation const &t;
  Polyhedron_rt const iso_value;
  T2P_converter const &t2p_converter;
  P2T_converter const &p2t_converter;
  Triang_vertex_map triang_vertex_signs;
};

CGAL_END_NAMESPACE

#endif // SKIN_SURFACE_REFINEMENT_TRAITS_H
