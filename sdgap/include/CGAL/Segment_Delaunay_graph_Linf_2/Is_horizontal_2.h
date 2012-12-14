
#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_IS_HORIZONTAL_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_IS_HORIZONTAL_2_H

#include <CGAL/Segment_Delaunay_graph_Linf_2/basic.h>

namespace CGAL {

namespace SegmentDelaunayGraphLinf_2 {

//-----------------------------------------------------------------------
//                           Is horizontal
//-----------------------------------------------------------------------

template< class K >
class Is_horizontal_2
{
private:
  typedef typename K::Site_2              Site_2;
  typedef typename K::Segment_2           Segment_2;
  typedef typename K::Compare_y_2         Compare_y_2;
  typedef typename K::Boolean             Boolean;
  
  Compare_y_2 compare_y_2;

public:

  Boolean operator()(const Site_2& s) const
  {
    CGAL_precondition( s.is_segment() );
    return compare_y_2(s.source(),s.target()) == EQUAL;
  }

};

} //namespace SegmentDelaunayGraphLinf_2

} //namespace CGAL


#endif //  CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_IS_HORIZONTAL_2_H
