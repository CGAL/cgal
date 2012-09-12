
namespace CGAL {

/*!
  \ingroup PkgArrangement2

  The traits class `Arr_rational_function_traits_2` is a model of the `ArrangementTraits_2` 
  concept. It handles bounded and unbounded arcs of rational functions, 
  referred to as <I>rational arcs</I> (in particular, such an arc may 
  correspond to the entire graph of a rational function), and enables the 
  construction and maintenance of arrangements of such arcs. 

  A rational function \f$ y = \frac{P(x)}{Q(x)}\f$ is defined by two polynomials 
  \f$ P\f$ and \f$ Q\f$ of arbitrary degrees. 
  If \f$ Q(x) = 1\f$ then the function is a simple polynomial function. 
  Usually the domain is \f$ \R\f$ but the function may also be 
  restricted to a bounded interval \f$ [x_{\rm min}, x_{\rm max}]\f$ 
  or defined over a ray \f$ (-\infty, x_{\rm max}]\f$ or over \f$ [x_{\rm min}, \infty)\f$. 
  Rational functions are represented by the nested type `Curve_2`. 
  Note that a rational function may be not continuous since roots of \f$ Q\f$ induce 
  vertical asymptotes, which would contradict the notion of an \f$ x\f$-monotone curve 
  as it is introduced by the `ArrangementTraits_2` concept. 
  Thus, continuous portions of rational functions are represented by the nested 
  type `X_monotone_curve_2`, which is different from `Curve_2`. 
  Constructors for both classes are provided by the traits. 
  A `Curve_2` may be split up into several `X_monotone_curve_2` 
  using `Make_x_monotone_2`. 

  The template parameter of the traits must be a model of the 
  concept `AlgebraicKernel_d_1`. 
  A rational function is then represented by two polynomials \f$ P\f$ and \f$ Q\f$ of type 
  `AlgebraicKernel_d_1::Polynomial_1`. 
  A point is represented by a rational function and its \f$ x\f$-coordinate, which is 
  of type `AlgebraicKernel_d_1::Algebraic_real_1`. 
  Note that an explicit representation of the \f$ y\f$-coordinate is only computed upon 
  request, which can be a rather costly operation. 

  The constructed rational functions are cached by the traits class. 
  The cache is local to each traits class object. 
  It is therefore necessary to construct the curves using the constructor 
  objects provided by member functions of the traits class. 
  Moreover, a curve must only be used with its own traits. 
  The cache is automatically cleaned up from time to time. 
  The amortized clean up costs are constant. However, there is also a 
  separate member function that cleans up the cache on demand. 

  While `Arr_rational_function_traits_2` models the concept 
  `ArrangementDirectionalXMonotoneTraits_2`, the implementation of 
  the `Arr_mergeable_2` operation does not enforce the input curves 
  to have the same direction as a precondition. Moreover, `Arr_rational_function_traits_2` 
  supports the merging of curves of opposite directions. 

  \models ::ArrangementTraits_2 
  \models ::ArrangementDirectionalXMonotoneTraits_2 
  \models ::ArrangementOpenBoundaryTraits_2 

  CONVERROR 5 nested classes missing 

*/
template< typename AlgebraicKernel_d_1 >
class Arr_rational_function_traits_2 {
public:
  
  /// \name Types 
  /// @{

  /*! 

   */ 
  typedef AlgebraicKernel_d_1 Algebraic_kernel_d_1; 

  /*! 

   */ 
  typedef AlgebraicKernel_d_1::Coefficient Coefficient; 

  /*! 

   */ 
  typedef AlgebraicKernel_d_1::Polynomial_1 Polynomial_1; 

  /*! 

   */ 
  typedef AlgebraicKernel_d_1::Algebraic_real_1 Algebraic_real_1; 

  /*! 

   */ 
  typedef AlgebraicKernel_d_1::Bound Bound; 

  /// @} 

  /// \name Creation 
  /// @{

  /*! 
    constructs an empty traits that uses the kernel pointed by `kernel` 
    for performing algebraic operations. 
  */ 
  Arr_rational_function_traits_2<AlgebraicKernel_d_1>(const Algebraic_kernel_d_1* kernel); 

  /// @} 

  /// \name Operations 
  /// @{

  /*! 
    Returns an instance of `Construct_curve_2`. 
  */ 
  Construct_curve_2 construct_curve_2_object() const; 

  /*! 
    Returns an instance of `Construct_x_monotone_curve_2`. 
  */ 
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object() const; 

  /*! 
    Deletes all curves from the cache that exist only there. 
  */ 
  void cleanup_cache() const; 

  /*! 
    Returns a pointer to the used algerbaic kernel object. 
  */ 
  const Algebraic_kernel_d_1* algebraic_kernel_d_1() const; 

  /// @}

}; /* end Arr_rational_function_traits_2 */
} /* end namespace CGAL */
