#ifndef POLYGON_TRIANGULATION_TRAITS_2
#define POLYGON_TRIANGULATION_TRAITS_2

#include <CGAL/Pm_segment_traits_2.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/Monotone_polygon_triangulation_traits_2.h>

#include <CGAL/Circulator_project.h>

namespace CGAL {

template < class Node, class Result>
struct Project_halfedge_source_point {
  typedef Node                   argument_type;
  typedef Result                 result_type;
  const Result& operator()( const Node& x)  { return x.source()->point(); }
  Result&       operator()( Node& x)        { return x.source()->point(); }
};

template <class Ccb_halfedge_circulator>
class Ccb_halfedge_bidirectional_circulator : public Ccb_halfedge_circulator
{
  typedef Ccb_halfedge_circulator Base;
  typedef Ccb_halfedge_bidirectional_circulator<Ccb_halfedge_circulator> Self;

 public:
  Ccb_halfedge_bidirectional_circulator() : Base() {}
  Ccb_halfedge_bidirectional_circulator(const Base& fc) : Base(fc) {}

  Self& operator--() {
    Base prev;
    //if( (*this)->source()->degree() <= 2)
    //  prev = (*this)->twin()->next_halfedge()->twin();
    //else {
      prev = (*this)->twin();
      while( prev->next_halfedge() != *this)
        prev = prev->next_halfedge()->twin();
    //}
    *this = prev;
    return *this;
  }

  Self operator--(int) {
    Self tmp = *this;
    --*this;
    return tmp;
  }
};


template <class Kernel_>
class Polygon_triangulation_traits_2
{
 private:
  typedef Kernel_ Kernel;

 public:
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Triangle_2 Triangle_2;

  typedef Pm_segment_traits_2<Kernel> Pm_segment_traits_2;
  typedef Pm_default_dcel<Pm_segment_traits_2> Dcel;
  typedef Planar_map_2< Dcel, Pm_segment_traits_2> Planar_map_2;

  typedef Circulator_project<
    Ccb_halfedge_bidirectional_circulator<
      typename Planar_map_2::Ccb_halfedge_circulator>,
    Project_halfedge_source_point<
      typename Planar_map_2::Halfedge, Point_2>,
    Point_2&,
    Point_2*
  > Circulator_project;
  typedef std::pair< Circulator_project, Circulator_project> Diagonal;

  typedef Partition_traits_2<Diagonal, Kernel> Partition_traits_2;

  typedef Monotone_polygon_triangulation_traits_2<Diagonal, Kernel>
    Monotone_polygon_triangulation_traits_2;
};

}

#endif // POLYGON_TRIANGULATION_TRAITS_2
