/*!
\ingroup PkgConeBasedSpannersConcepts
\cgalConcept

The functors provided in this package for constructing cones and cone-based spanners
all have a template parameter `Traits`. To specify the requirements
for the models (typically the kernels from CGAL) that can be passed to 
this parameter, we document a concept called `ConeBasedSpannerTraits_2` here. 
Basically, this concept specifies all the types and primitives (predicates and construction objects) 
that the model should have to make the functors work properly.
It is recommended that if you want to construct the cones or the cone-based spanners
exactly, you should use the kernel `CGAL::Exact_predicates_exact_constructions_kernel_with_root_of`;
and if you want to construct them inexactly, you should use the kernel
`CGAL::Exact_predicates_inexact_constructions`.

\cgalHasModel `CGAL::Exact_predicates_exact_constructions_kernel_with_root_of`
\cgalHasModel `CGAL::Exact_predicates_inexact_constructions`

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
  by the iterator range. It is needed in calculating cone boundaries exactly.
*/ 
NT CGAL::root_of(int k, InputIterator begin, InputIterator end); 

/*!
  This function should return the square root of the argument `x`. 
*/ 
NT CGAL::sqrt(const NT &  x); 

/*!
  This function should return the bisector of the two lines l1 and l2. 
  And the bisector should have the direction of the vector which is the sum of 
  the normalized directions of the two lines. 
  It is needed in constructing Theta graphs, not in constructing Yao graphs.
*/ 
CGAL::Line_2<Kernel> CGAL::bisector(const CGAL::Line_2<Kernel> &  l1,  
		  const CGAL::Line_2<Kernel> &  l2);

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
