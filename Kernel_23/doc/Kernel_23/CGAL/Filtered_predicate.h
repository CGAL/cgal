
namespace CGAL {

/*!
\ingroup kernel_classes

`Filtered_predicate` is an adaptor for predicate function objects that allows 
one to produce efficient and exact predicates. It is used to build 
`CGAL::Filtered_kernel<CK>` and can be used for other 
predicates too. 

`EP` is the exact but supposedly slow predicate that is able to evaluate 
the predicate correctly. It will be called only when the filtering predicate, 
`FP`, cannot compute the correct result. This failure of `FP` must be 
done by throwing an exception. 

To convert the geometric objects that are the arguments of the predicate, 
we use the function objects `C2E` and `C2F`, which must be of the form 
`Cartesian_converter` or `Homogeneous_converter`. 

\cgalHeading{Example}

The following example defines an efficient and exact version of the 
orientation predicate over three points using the Cartesian representation 
with double coordinates and without reference counting 
(`Simple_cartesian::Point_2`). 
Of course, the orientation predicate can already be found in the kernel, but 
you can follow this example to filter your own predicates. 
It uses the fast but inexact predicate based on interval arithmetic for 
filtering and the slow but exact predicate based on multi-precision floats 
when the filtering predicate fails. 

\cgalExample{Filtered_kernel/Filtered_predicate.cpp} 

*/
template< typename EP, typename FP, typename C2E, typename C2F >
class Filtered_predicate {
public:

/// \name Types 
/// @{

/*!
The return type of the function operators. 
It must also be the same type as `EP::result_type`. 
*/ 
typedef FP::result_type result_type; 

/// @} 

/// \name Creation 
/// @{

/*!
%Default constructor. 
*/ 
Filtered_predicate<EP, FP, C2E, C2F>(); 

/// @} 

/// \name Operations 
/// Similar function operators are defined for up to 7 arguments.
/// @{

/*!
The unary function operator for unary predicates. 
*/ 
template <class A1> result_type operator()(A1 a1); 

/*!
The binary function operator for binary predicates. 
*/ 
template <class A1, class A2> 
result_type operator()(A1 a1, A2 a2); 

/// @}

}; /* end Filtered_predicate */
} /* end namespace CGAL */
