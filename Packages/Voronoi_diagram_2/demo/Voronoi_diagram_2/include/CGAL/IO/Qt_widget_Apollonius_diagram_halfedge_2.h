#ifndef CGAL_QT_WIDGET_APOLLONIUS_DIAGRAM_HALFEDGE_2_H
#define CGAL_QT_WIDGET_APOLLONIUS_DIAGRAM_HALFEDGE_2_H 1

#ifdef CGAL_USE_QT

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/Apollonius_graph_constructions_C2.h>
#include <CGAL/Hyperbola_segment_2.h>
#include <CGAL/Hyperbola_ray_2.h>
#include <CGAL/Hyperbola_2.h>


CGAL_BEGIN_NAMESPACE

template<class VDA>
class Apollonius_diagram_halfedge_2
  : public VDA::Halfedge
{
 protected:
  typedef VDA                                        Voronoi_diagram;
  typedef typename Voronoi_diagram::Delaunay_graph   Apollonius_graph_2;
  typedef typename Voronoi_diagram::Halfedge         Base;
  typedef typename Base::Delaunay_edge               Delaunay_edge;

  typedef typename Voronoi_diagram::Voronoi_traits::Site_2  Site_2;

 public:
  Apollonius_diagram_halfedge_2() : Base() {}
  Apollonius_diagram_halfedge_2(const Base& e)
    : Base(e), is_conflict(false) {}
  Apollonius_diagram_halfedge_2(const Delaunay_edge& e, bool inf,
				const Site_2& s)
    : Base(), is_conflict(true), e_(e), inf_(inf), s_(s) {}

  void draw(Qt_widget& qt_w) const
  {
    if ( is_conflict ) { return; }

    typedef typename Apollonius_graph_2::Geom_traits     Geom_traits;
    typedef typename Geom_traits::Assign_2               Assign_2;
    typedef typename Geom_traits::Segment_2              Segment_2;
    typedef typename Geom_traits::Ray_2                  Ray_2;
    typedef typename Geom_traits::Line_2                 Line_2;
    typedef Hyperbola_segment_2<Geom_traits>             Hyperbola_segment_2;
    typedef Hyperbola_ray_2<Geom_traits>                 Hyperbola_ray_2;
    typedef Hyperbola_2<Geom_traits>                     Hyperbola_2;

    Assign_2 assign = Geom_traits().assign_2_object();
    Hyperbola_segment_2 hs;
    Hyperbola_ray_2 hr;
    Hyperbola_2 h;
    Segment_2 s;
    Ray_2 r;
    Line_2 l;

    Object o;
    if ( this->has_source() && this->has_target() ) {
      Construct_Apollonius_bisector_segment_2<Geom_traits> c_seg;
      o = c_seg(this->down()->site(),
		this->up()->site(),
		this->left()->site(),
		this->right()->site());
    } else if ( this->has_source() && !this->has_target() ) {
      Construct_Apollonius_bisector_ray_2<Geom_traits> c_ray;
      o = c_ray(this->down()->site(),
		this->up()->site(),
		this->left()->site());
    } else if ( !this->has_source() && this->has_target() ) {
      Construct_Apollonius_bisector_ray_2<Geom_traits> c_ray;
      o = c_ray(this->up()->site(),
		this->down()->site(),
		this->right()->site());
    } else {
      CGAL_assertion( !this->has_source() && !this->has_target() );
      Construct_Apollonius_bisector_2<Geom_traits> c_bis;
      o = c_bis(this->up()->site(),
		this->down()->site());
    }

    // fix this and use the output operators...
    if      ( assign(hs,o) )   hs.draw(qt_w);
    else if ( assign(hr,o) )   hr.draw(qt_w);
    else if ( assign(h, o) )   h.draw(qt_w);
    else if ( assign(s, o) )   qt_w << s;
    else if ( assign(r, o) )   qt_w << r;
    else if ( assign(l, o) )   qt_w << l;
  }

private:
  bool is_conflict;
  Delaunay_edge e_;
  bool inf_;
  Site_2 s_;
};

template<class VDA>
Qt_widget& operator<<(Qt_widget& qt_w,
		      const Apollonius_diagram_halfedge_2<VDA>& e)
{
  e.draw(qt_w);
  return qt_w;
}


CGAL_END_NAMESPACE

#endif // CGAL_USE_QT

#endif // CGAL_QT_WIDGET_APOLLONIUS_DIAGRAM_HALFEDGE_2_H
