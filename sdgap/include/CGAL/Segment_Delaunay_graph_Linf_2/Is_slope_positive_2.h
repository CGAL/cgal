
#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_IS_SLOPE_POSITIVE_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_IS_SLOPE_POSITIVE_2_H

#include <CGAL/Segment_Delaunay_graph_Linf_2/basic.h>

namespace CGAL {

namespace SegmentDelaunayGraphLinf_2 {

//-----------------------------------------------------------------------
//                           Is slope positive
//-----------------------------------------------------------------------

template< class K >
class Is_slope_positive_2
{
private:
  typedef typename K::Site_2              Site_2;
  typedef typename K::Segment_2           Segment_2;
  typedef typename K::Is_horizontal_2     Is_horizontal_2;
  typedef typename K::Is_vertical_2       Is_vertical_2;
//  typedef typename K::Compare_x_2         Compare_x_2;
//  typedef typename K::Compare_y_2         Compare_y_2;
  typedef typename K::Boolean             Boolean;
  
//  Compare_x_2 compare_x_2;
//  Compare_y_2 compare_y_2;
  Is_horizontal_2 is_horizontal_2;
  Is_vertical_2 is_vertical_2;

public:

  Boolean operator()(const Site_2& s) const
  {
    CGAL_precondition( s.is_segment() );
    CGAL_precondition( !is_horizontal_2(s) && !is_vertical_2(s) );
    return CGAL::sign(s.target().x() - s.source().x()) ==
           CGAL::sign(s.target().y() - s.source().y());
  }

};

} //namespace SegmentDelaunayGraphLinf_2

} //namespace CGAL


#endif //  CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_IS_SLOPE_POSITIVE_2_H
