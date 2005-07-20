#ifndef CGAL_QT_WIDGET_VORONOI_DIAGRAM_HALFEDGE_2_H
#define CGAL_QT_WIDGET_VORONOI_DIAGRAM_HALFEDGE_2_H 1

#ifdef CGAL_USE_QT

#include <CGAL/IO/Qt_widget.h>

CGAL_BEGIN_NAMESPACE

template<class VDA>
class Voronoi_diagram_halfedge_2
  : public VDA::Halfedge
{
 protected:
  typedef VDA                                        Voronoi_diagram;
  typedef typename Voronoi_diagram::Delaunay_graph   Delaunay_triangulation_2;
  typedef typename Voronoi_diagram::Halfedge         Base;

 public:
  Voronoi_diagram_halfedge_2() : Base() {}
  Voronoi_diagram_halfedge_2(const Base& e) : Base(e) {}

  void draw(Qt_widget& qt_w) const
  {
    typedef typename Delaunay_triangulation_2::Geom_traits   Geom_traits;
    typedef typename Geom_traits::Point_2                    Point_2;
    typedef typename Geom_traits::Segment_2                  Segment_2;
    typedef typename Geom_traits::Line_2                     Line_2;

    if ( this->has_source() && this->has_target() ) {
      typename Geom_traits::Construct_circumcenter_2 circumcenter;
      Point_2 c1 = circumcenter(this->down()->point(),
				this->up()->point(),
				this->left()->point());
      Point_2 c2 = circumcenter(this->up()->point(),
				this->down()->point(),
				this->right()->point());
      qt_w << Segment_2(c1, c2);
    } else if ( this->has_source() && !this->has_target() ) {
      typename Geom_traits::Construct_circumcenter_2 circumcenter;
      typename Geom_traits::Construct_bisector_2     c_bis;
      typename Geom_traits::Construct_ray_2          c_ray;
      Point_2 c = circumcenter(this->down()->point(),
			       this->up()->point(),
			       this->left()->point());
      Line_2 l = c_bis(this->up()->point(), this->down()->point());
      qt_w << c_ray(c, l);
    } else if ( !this->has_source() && this->has_target() ) {
      typename Geom_traits::Construct_circumcenter_2 circumcenter;
      typename Geom_traits::Construct_bisector_2     c_bis;
      typename Geom_traits::Construct_ray_2          c_ray;
      Point_2 c = circumcenter(this->up()->point(),
			       this->down()->point(),
			       this->right()->point());
      Line_2 l = c_bis(this->down()->point(), this->up()->point());
      qt_w << c_ray(c, l);
    } else {
      CGAL_assertion( !this->has_source() && !this->has_target() );
      typename Geom_traits::Construct_bisector_2 c_bis;
      qt_w << c_bis(this->up()->point(),
		    this->down()->point());
    }
  }
};

template<class VDA>
Qt_widget& operator<<(Qt_widget& qt_w,
		      const Voronoi_diagram_halfedge_2<VDA>& e)
{
  e.draw(qt_w);
  return qt_w;
}


CGAL_END_NAMESPACE

#endif // CGAL_USE_QT

#endif // CGAL_QT_WIDGET_VORONOI_DIAGRAM_HALFEDGE_2_H
