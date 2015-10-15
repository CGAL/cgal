/*!
\ingroup PkgConeBasedSpannersConcepts
\cgalConcept

The functors provided in this package for constructing
cone-based spanners are parameterized with a traits class `Traits`, which defines the 
primitives (predicates and construction objects) required by the functors.

\cgalHasModel any model of a \cgal %kernel.

*/

class ConeBasedSpannerTraits {
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
The affine transformation type. 
*/ 
typedef unspecified_type Aff_transformation_2; 

/*!
The polynomial type. When the cone angle \f$ \theta \f$ is in the form of 
\f$ 2\pi / n \f$, where \f$ n \f$ is a positive integer, \f$ \sin(\theta) \f$ 
and \f$ \cos(\theta) \f$ can be represented exactly by roots of polynomials.
*/ 
typedef unspecified_type Polynomial; 

/// @} 


/// \name Functions 
/// The following functions are needed to construct cone-based spanners.
/// @{

/*!
  This function should return the k-th real root of an univariate polynomial, which is defined 
  by the iterator range. It is needed in calculating cone boundaries exactly.
  It is not needed if the cone bounaries are calculated inexactly, which uses
  sin() and cos() instead.
*/ 
NT CGAL::root_of(int k, InputIterator begin, InputIterator end); 

/*!
  This function should return the square root of the argument `x`. 
  It is defined if the argument type is a model of the `FieldWithSqrt` concept.
  It is needed in calculating cone boundaries exactly.
*/ 
NT CGAL::sqrt(const NT &  x); 

/*!
  This function should compute the sine value of arg (measured in radians). 
  It is needed in calculating cone boundaries inexactly.
*/ 
long double sin( long double arg );

/*!
  This function should compute the cosine value of arg (measured in radians). 
  It is needed in calculating cone boundaries inexactly.
*/ 
long double cos( long double arg );

/*!
  This function should return the bisector of the two lines l1 and l2. 
  And the bisector should have the direction of the vector which is the sum of 
  the normalized directions of the two lines. 
  This function requires that Kernel::RT supports the sqrt() operation. 
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

}; /* end ConeBasedSpannerTraits */
