#ifndef CGAL_CURVED_KERNEL_GLOBAL_FUNCTIONS_ON_LINE_2_H
#define CGAL_CURVED_KERNEL_GLOBAL_FUNCTIONS_ON_LINE_2_H

namespace CGAL {

template< class CK >
inline
typename CK::Polynomial_1_2
get_equation(const typename CK::Line_2 & l)
{
  return CK().get_equation_object()(l);
}

template< class CK >
inline
typename CK::Line_2
construct_line_2(const typename CK::Polynomial_1_2 & eq)
{
  return CK().construct_line_2_object()(eq);
}

template< class CK, class OutputIterator>
inline
OutputIterator
construct_intersections_2( const typename CK::Line_2 & l,
			   const typename CK::Circle_2 & c,
			   OutputIterator res )
{
  return CK().construct_intersections_2_object()(l,c,res);
}


} // namespace CGAL
#endif // CGAL_CURVED_KERNEL_GLOBAL_FUNCTIONS_ON_LINE_2_H
