#ifndef CGAL_LINE_ARC_2_H
#define CGAL_LINE_ARC_2_H

namespace CGAL {

template <class CurvedKernel> 
class Line_arc_2 
  : public CurvedKernel::Kernel_base::Line_arc_2
{
  typedef typename CurvedKernel::FT                        FT;
  typedef typename CurvedKernel::RT                        RT;
  typedef typename CurvedKernel::Linear_kernel::Point_2    Point_2;
  typedef typename CurvedKernel::Line_2                    Line_2;
  typedef typename CurvedKernel::Circle_2                  Circle_2;
  typedef typename CurvedKernel::Circular_arc_point_2   Circular_arc_point_2;
  typedef typename CurvedKernel::Segment_2                 Segment_2;

  typedef typename CurvedKernel::Kernel_base::Line_arc_2 RLine_arc_2;
public:
  typedef  RLine_arc_2 Rep;
  typedef  CurvedKernel   R; 

 const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

   Line_arc_2()
     : RLine_arc_2(typename R::Construct_line_arc_2()())
   {}

   Line_arc_2(const Line_2 &support,
	      const Circle_2 &c1,const bool b1,
	      const Circle_2 &c2,const bool b2)
     : RLine_arc_2(typename R::Construct_line_arc_2()(support, c1, b1, c2, b2))
   {}

   Line_arc_2(const Line_2 &support,
	       const Line_2 &l1,
	       const Line_2 &l2)
     : RLine_arc_2(typename R::Construct_line_arc_2()(support, l1, l2))
   {}

   Line_arc_2(const Line_2 &support,
	       const Circular_arc_point_2 &p1,
	       const Circular_arc_point_2 &p2)
     : RLine_arc_2(typename R::Construct_line_arc_2()(support, p1, p2))
   {}

   Line_arc_2(const Segment_2 &s)
     : RLine_arc_2(typename R::Construct_line_arc_2()(s))
   {}
   
    Line_arc_2(const Point_2 &p1,
	       const Point_2 &p2)
      : RLine_arc_2(typename R::Construct_line_arc_2()(p1, p2))
   {}

    Line_arc_2(const RLine_arc_2 &a )
     : RLine_arc_2(a)
   {}

  const Circular_arc_point_2 & source() const
  {
        return typename R::Construct_source_vertex_2()(*this);
  }

  const Circular_arc_point_2 & target() const
  {
        return typename R::Construct_target_vertex_2()(*this);
  }

  const Circular_arc_point_2 & left() const
  {
        return typename R::Construct_min_vertex_2()(*this);
  }

  const Circular_arc_point_2 & right() const
  {
        return typename R::Construct_max_vertex_2()(*this);
  }


  const Line_2 & supporting_line() const
  {
        return typename R::Construct_supporting_line_2()(*this);
  }
  

  bool is_vertical() const
  {
      return typename R::Is_vertical_2()(*this);
  }
    
  Bbox_2  bbox() const
  {
        return typename R::Construct_bbox_2()(*this);
  }

 };

template < typename CurvedKernel >
inline
bool
operator==(const Line_arc_2<CurvedKernel> &p,
           const Line_arc_2<CurvedKernel> &q)
{
  return CurvedKernel().equal_2_object()(p, q);
}

template < typename CurvedKernel >
inline
bool
operator!=(const Line_arc_2<CurvedKernel> &p,
           const Line_arc_2<CurvedKernel> &q)
{
  return ! (p == q);
}










} // namespace CGAL

#endif // CGAL_LINE_ARC_2_H

   
