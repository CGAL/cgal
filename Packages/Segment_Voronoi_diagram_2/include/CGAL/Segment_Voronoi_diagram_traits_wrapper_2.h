#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_TRAITS_WRAPPER_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_TRAITS_WRAPPER_2_H

CGAL_BEGIN_NAMESPACE


template<class Gt_base>
class Segment_Voronoi_diagram_traits_wrapper_2 : public Gt_base
{
private:
  typedef Segment_Voronoi_diagram_traits_wrapper_2<Gt_base> Self;

public:
  struct Triangle_2 {};

  Segment_Voronoi_diagram_traits_wrapper_2() {}

  Segment_Voronoi_diagram_traits_wrapper_2
  (const Self&) {}

  Segment_Voronoi_diagram_traits_wrapper_2
  operator=(const Self&) {
    return (*this);
  }

  Segment_Voronoi_diagram_traits_wrapper_2(const Gt_base&) {}
};




CGAL_END_NAMESPACE


#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_TRAITS_WRAPPER_2_H
