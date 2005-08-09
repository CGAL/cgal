
// Author : Constantinos Tsirogiannis

#ifndef CGAL_CURVED_KERNEL_CONVERTER_H
#define CGAL_CURVED_KERNEL_CONVERTER_H 

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Algebraic_kernel_converter.h>

// TODO :
// - we should have a better default than Cartesian_converter.

CGAL_BEGIN_NAMESPACE

template < class C1, class C2,
	   class LK_converter = Cartesian_converter<typename C1::Linear_kernel,
                                                    typename C2::Linear_kernel>,
	   class AK_converter = Algebraic_kernel_converter<typename C1::Algebraic_kernel,
							 typename C2::Algebraic_kernel > >
class Curved_kernel_converter 
{
public:

	typedef C1                       		      Source_kernel;
	typedef C2                          		      Target_kernel;
	typedef typename C1::Linear_kernel  		      L1;
	typedef typename C2::Linear_kernel  		      L2;
	typedef LK_converter                		      Linear_kernel_converter;
	typedef AK_converter                		      Algebraic_kernel_converter;
	typedef typename Linear_kernel_converter::Number_type_converter       RT_type_converter;
	typedef typename Algebraic_kernel_converter::Root_of_type_converter   Root_of_type_converter;

	typename C2::Circle_2
	operator()(const typename C1::Circle_2 &a) const
	{
		return Linear_kernel_converter()(a);
	}

	typename C2::Circular_arc_endpoint_2
	operator()(const typename C1::Circular_arc_endpoint_2 &a) const
	{
		return typename C2::Circular_arc_endpoint_2( typename C2::Circular_arc_endpoint_2::Numeric_point_2(
							     Root_of_type_converter()( a.x() ),
							     Root_of_type_converter()( a.y() )
	   	  )
	   );
	}

	typename C2::Circular_arc_2
	operator()(const typename C1::Circular_arc_2 &a) const
	{
		return typename C2::Circular_arc_2 ( operator()( a.supporting_circle() ),
		                     		     operator()( a.source() ),
		                     		     operator()( a.target() ) );
	}
}; 

CGAL_END_NAMESPACE

#endif // CGAL_CURVED_KERNEL_CONVERTER_H
