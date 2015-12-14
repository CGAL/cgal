/*!
\ingroup PkgConeBasedSpannersConcepts
\cgalConcept

The concept `ConeBasedSpannerTraits_2` specifies the requirements
for the CGAL kernels that can be used in this package. 
Basically, this concept describes all the types and primitives (predicates and construction objects) 
that the CGAL kernel should have to make this package work properly.
It is recommended that if the cone-based spanners are to be constructed
exactly, the kernel `CGAL::Exact_predicates_exact_constructions_kernel_with_root_of`
should be used; and if they are to be constructed inexactly, the kernel
`CGAL::Exact_predicates_inexact_constructions_kernel` should be used.

\cgalHasModel `CGAL::Exact_predicates_exact_constructions_kernel_with_root_of`
\cgalHasModel `CGAL::Exact_predicates_inexact_constructions_kernel`

*/

class ConeBasedSpannerTraits_2 {
public:

/// \name Types 
/// @{

/*!
The point type. 
*/ 
typedef unspecified_type Point_2; 

/*!
The line type. 
*/ 
typedef unspecified_type Line_2; 

/*!
The direction type. 
*/ 
typedef unspecified_type Direction_2; 

/*!
The affine transformation type. This is needed to rotate directions.
*/ 
typedef unspecified_type Aff_transformation_2; 

/*!
The polynomial type. When the cone angle \f$ \theta \f$ is in the form of 
\f$ 2\pi / n \f$, where \f$ n \f$ is a positive integer, \f$ \sin(\theta) \f$ 
and \f$ \cos(\theta) \f$ can be represented exactly by roots of polynomials.
Thus, this polynomial type is needed to avoid the computation by calling
sin() and cos().
*/ 
typedef unspecified_type Polynomial; 

/// @} 


/// \name Functions 
/// The following functions are needed to construct cone-based spanners.
/// @{

/*!
  This function should return the k-th real root of an univariate polynomial, which is defined 
  by the iterator range, where 'begin' refers to the constant term. 
  It is needed in calculating cone boundaries exactly.
*/ 
NT CGAL::root_of(int k, InputIterator begin, InputIterator end); 

/*
  This function should return the square root of the argument `x`. 
  It is needed in calculating cone boundaries exactly.

NT CGAL::sqrt(const NT &  x); 
*/ 

/*!
  This functor should return the bisector of the two lines l1 and l2. 
  And the bisector should have the direction of the vector which is the sum of 
  the normalized directions of the two lines. 
  It is needed in constructing Theta graphs, not in constructing Yao graphs.
*/ 
Kernel::Line_2 Kernel::ConstructBisector_2::operator()  ( const Kernel::Line_2 &  l1,  
		  const Kernel::Line_2 &  l2 ); 

/*!
  This function should return CGAL::LARGER iff the signed distance of p and l 
  is larger than the signed distance of q and l, CGAL::SMALLER, iff it is smaller, 
  and CGAL::EQUAL iff both are equal. 
  It is needed in sorting the points on the plane based on a certain direction.
*/ 
Comparison_result  CGAL::compare_signed_distance_to_line(
		               const CGAL::Line_2<Kernel> &l, 
					   const CGAL::Point_2<Kernel> &p, 
					   const CGAL::Point_2<Kernel> &q);
/// @}

}; /* end ConeBasedSpannerTraits_2 */
