#ifndef CGAL_VORONOI_DIAGRAM_2_VORONOI_VERTEX_BASE_2_H
#define CGAL_VORONOI_DIAGRAM_2_VORONOI_VERTEX_BASE_2_H 1

#include <CGAL/Voronoi_diagram_adaptor_2/basic.h>

CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

template<class DG, class P, class S, class Voronoi_vertex>
class Voronoi_vertex_base_2
{
 private:
  typedef Voronoi_vertex_base_2<DG,P,S,Voronoi_vertex>   Self;

 public:
  typedef DG                                   Dual_graph;
  typedef typename Dual_graph::Geom_traits     Geom_traits;
  //  typedef typename Dual_graph::Vertex_handle   Vertex_handle;
  typedef P                                    Point_2;
  typedef S                                    Site_2;
  typedef Voronoi_vertex                       Voronoi_vertex_2;


  const Site_2& first() const { return s_[0]; }
  const Site_2& second() const { return s_[1]; }
  const Site_2& third() const { return s_[2]; }

  const Site_2& site(unsigned int i) const {
    CGAL_precondition( i <= 2 );
    return s_[i];
  }

#if 0
  bool operator==(const Voronoi_vertex_2& o) const {
    bool equal_found = false;
    unsigned int i_eq = 3;
    for (unsigned int i = 0; i < 3; i++) {
      if ( Geom_traits().equal_2_object()(s_[i], o.s_[0]) ) {
	equal_found = true;
	i_eq = i;
	break;
      }
    }

    if ( equal_found ) {
      if ( Geom_traits().equal_2_object()(s_[(i_eq+1)%3], o.s_[1]) &&
	   Geom_traits().equal_2_object()(s_[(i_eq+2)%3], o.s_[2]) )
      return true;
    }

    return false;
#if 0
    Point_2 p_this(*this);
    Point_2 p_other(o);
    // THIS MAY NOT BE EXACT
    return geom_traits().equal_2_object()(p_this, p_other);
#endif
  }

  bool operator!=(const Voronoi_vertex_2& o) const {
    return !(*this == o);
  }
#endif

 protected:
  void set_sites(const Site_2& s0, const Site_2& s1, const Site_2& s2) {
    s_[0] = s0;
    s_[1] = s1;
    s_[2] = s2;
  }

 protected:
  Site_2 s_[3];
};


CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE

#endif // CGAL_VORONOI_DIAGRAM_2_VORONOI_VERTEX_BASE_2_H
