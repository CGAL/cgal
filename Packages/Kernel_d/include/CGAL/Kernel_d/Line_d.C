#ifndef CGAL_LINE_D_C
#define CGAL_LINE_D_C
CGAL_BEGIN_NAMESPACE

template <class R> 
Line_d<R> Segment_d<R>::supporting_line() const
{ CGAL_assertion_msg((!is_degenerate()), 
  "Segment_d::supporting_line(): degenerate segment cannot be converted.");
  return Line_d<R>(Base(*this)); 
} 

template <class R>
Line_d<R> Ray_d<R>::supporting_line() const
{ return Line_d<R>(Base(*this)); } 

CGAL_END_NAMESPACE
#endif //CGAL_LINE_D_C

