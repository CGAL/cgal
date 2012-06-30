
#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_ORIENTATION_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_ORIENTATION_C2_H

#include <CGAL/Segment_Delaunay_graph_Linf_2/Basic_predicates_C2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Orientation_Linf_C2.h>

namespace CGAL {

namespace SegmentDelaunayGraphLinf_2 {

//-----------------------------------------------------------------------------



template<class K>
class Orientation_C2
  : private Basic_predicates_C2<K>
{
private:
  typedef Basic_predicates_C2<K>              Base;
  
public:
  typedef typename Base::Orientation          Orientation;

private:
  typedef typename Base::Point_2              Point_2;
  typedef typename Base::Segment_2            Segment_2;
  typedef typename Base::Site_2               Site_2;
  typedef typename Base::FT                   FT;
  typedef typename Base::RT                   RT;

  typedef typename Base::Comparison_result    Comparison_result;
  typedef typename Base::Oriented_side        Oriented_side;
  typedef typename Base::Sign                 Sign;

  typedef Orientation_Linf_C2<K>      Orientation_Linf; 
  //typedef typename K::Kernel::Orientation_2   Orientation_2;
  

  typedef typename K::Intersections_tag       ITag;

  //-------------------------------------------------------------

  Orientation predicate(const Site_2& p, const Site_2& q,
			const Site_2& r, const Tag_false&) const
  {
    return Orientation_Linf()(p.point(), q.point(), r.point());
  }

  Orientation predicate(const Site_2& p, const Site_2& q,
			const Site_2& r, const Tag_true&) const
  {
    return predicate(p, q, r, Tag_false());
  }

public:
  typedef Orientation                  result_type;
  typedef Site_2                       argument_type;

  Orientation operator()(const Site_2& p, const Site_2& q,
			 const Site_2& r) const
  {
    CGAL_precondition( p.is_point() && q.is_point() && r.is_point() );

    return predicate(p, q, r, ITag());
  }
};


//-----------------------------------------------------------------------------

} //namespace SegmentDelaunayGraphLinf_2

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_ORIENTATION_C2_H
