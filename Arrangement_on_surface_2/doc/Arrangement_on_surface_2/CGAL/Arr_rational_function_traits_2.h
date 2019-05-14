
namespace CGAL {

/*!
  \ingroup PkgArrangement2TraitsClasses

  The traits class `Arr_rational_function_traits_2` is a model of the `ArrangementTraits_2` 
  concept. It handles bounded and unbounded arcs of rational functions, 
  referred to as <i>rational arcs</i> (in particular, such an arc may correspond to the
  entire graph of a rational function). It supports bounded and
  unbounded arcs. Thus, it is also a model of the concept
  `ArrangementOpenBoundaryTraits_2`. The traits class enables
  the construction and maintenance of arrangements of such arcs. 

  A rational function \f$ y = \frac{P(x)}{Q(x)}\f$ is defined by two polynomials 
  \f$ P\f$ and \f$ Q\f$ of arbitrary degrees. 
  If \f$ Q(x) = 1\f$ then the function is a simple polynomial function. 
  Usually the domain is \f$ \mathbb{R}\f$ but the function may also be 
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
  the `Are_mergeable_2` operation does not enforce the input curves 
  to have the same direction as a precondition. Moreover, `Arr_rational_function_traits_2` 
  supports the merging of curves of opposite directions. 

  \cgalModels `ArrangementTraits_2`
  \cgalModels `ArrangementDirectionalXMonotoneTraits_2`
  \cgalModels `ArrangementOpenBoundaryTraits_2`
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


/*!


Functor to construct a `Curve_2`. To enable caching the class is not 
default constructible and must be obtained via the function 
`construct_curve_2_object()`, which is a member of the traits. 

\cgalModels `Assignable`
\cgalModels `CopyConstructible`
\cgalModels `AdaptableBinaryFunction`
\cgalModels `AdaptableUnaryFunction`

*/
class Construct_curve_2 {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef AlgebraicKernel_d_1::Polynomial_1 Polynomial_1; 

/*!

*/ 
typedef AlgebraicKernel_d_1::Algebraic_real_1 Algebraic_real_1; 

/*!

*/ 
typedef Arr_rational_function_traits_2<AlgebraicKernel_d_1>::Curve_2 result_type; 

/*!

*/ 
typedef Polynomial_1 argument_type; 

/*!

*/ 
typedef Polynomial_1 first_argument_type; 

/*!

*/ 
typedef Polynomial_1 second_argument_type; 

/// @} 

/// \name Operations 
/// @{

/*!
Constructs a curve representing the polynomial function \f$ y = P(x)\f$. 
*/ 
Curve_2 operator()(Polynomial_1 P) const; 

/*!
Constructs a curve representing the polynomial function \f$ y = P(x)\f$. 
The function is defined over the interval \f$ [x,+\infty)\f$ if \f$ right\f$ is true 
and \f$ (-\infty,x]\f$ otherwise. 
*/ 
Curve_2 operator()(Polynomial_1 P, const Algebraic_real_1& x, bool right) const; 

/*!
Constructs a curve representing the polynomial function \f$ y = P(x)\f$. 
The function is defined over the interval \f$ [lower,upper]\f$. 
*/ 
Curve_2 operator()(Polynomial_1 P, const Algebraic_real_1& lower, const Algebraic_real_1& upper) const; 

/*!
Constructs a curve representing the rational function \f$ y = P(x)/Q(x)\f$. 
*/ 
Curve_2 operator()(Polynomial_1 P, Polynomial_1 Q) const; 

/*!
Constructs a curve representing the rational function \f$ y = P(x)/Q(x)\f$. 
The function is defined over the interval \f$ I=[x,+\infty)\f$ if \f$ right\f$ is 
true and \f$ I=(-\infty,x]\f$ otherwise. 
*/ 
Curve_2 operator()(Polynomial_1 P, Polynomial_1 Q, const Algebraic_real_1& x, bool right) const; 

/*!
Constructs a curve representing the rational function \f$ y = P(x)/Q(x)\f$. 
The function is defined over the interval \f$ I=[lower,upper]\f$. 
*/ 
Curve_2 operator()(Polynomial_1 P, Polynomial_1 Q, const Algebraic_real_1& lower, const Algebraic_real_1& upper) const; 

/*!
Constructs a curve representing the polynomial function \f$ y = P(x)\f$, where 
the coefficients of \f$ P\f$ are given in the range `[begin,end)`. 
*/ 
template <typename InputIterator> 
Curve_2 operator()(InputIterator begin, InputIterator end) const; 

/*!
Constructs a curve representing the polynomial function \f$ y = P(x)\f$, where 
the coefficients of \f$ P\f$ are given in the range `[begin,end)`. The 
function is defined over the interval \f$ [x,+\infty)\f$ if \f$ right\f$ is true 
and \f$ (-\infty,x]\f$ otherwise. 
*/ 
template <typename InputIterator> 
Curve_2 operator()(InputIterator begin, InputIterator end, 
const Algebraic_real_1& x, bool right) const; 

/*!
Constructs a curve representing the polynomial function \f$ y = P(x)\f$, where 
the coefficients of \f$ P\f$ are given in the range `[begin,end)`. The 
function is defined over the interval \f$ [lower,upper]\f$. 
*/ 
template <typename InputIterator> 
Curve_2 operator()(InputIterator begin, InputIterator end, 
const Algebraic_real_1& lower, 
const Algebraic_real_1& upper) const; 

/*!
Constructs a curve representing the rational function \f$ y = P(x)/Q(x)\f$, 
where the coefficients of \f$ P\f$ and \f$ Q\f$ are given in the ranges 
`[begin_numer,end_numer)` and `[begin_denom,end_denom)`, 
respectively. 
*/ 
template <typename InputIterator> 
Curve_2 operator()(InputIterator begin_numer, InputIterator end_numer, 
InputIterator begin_denom, InputIterator end_denom) const; 

/*!
Constructs a curve representing the rational function \f$ y = P(x)/Q(x)\f$, 
where the coefficients of \f$ P\f$ and \f$ Q\f$ are given in the ranges 
`[begin_numer,end_numer)` and `[begin_denom,end_denom)`, 
respectively. The function is defined over the interval \f$ I=[x,+\infty)\f$ 
if \f$ right\f$ is true and \f$ I=(-\infty,x]\f$ otherwise. 
*/ 
template <typename InputIterator> 
Curve_2 operator()(InputIterator begin_numer, InputIterator end_numer, 
InputIterator begin_denom, InputIterator end_denom, 
const Algebraic_real_1& x, bool right) const; 

/*!
Constructs a curve representing the rational function \f$ y = P(x)/Q(x)\f$, 
where the coefficients of \f$ P\f$ and \f$ Q\f$ are given in the ranges 
`[begin_numer,end_numer)` and `[begin_denom,end_denom)`, 
respectively. The function is defined over the interval \f$ I=[lower,upper]\f$. 
*/ 
template <typename InputIterator> 
Curve_2 operator()(InputIterator begin_numer, InputIterator end_numer, 
InputIterator begin_denom, InputIterator end_denom, 
const Algebraic_real_1& lower, 
const Algebraic_real_1& upper) const; 

/// @}

}; /* end Arr_rational_function_traits_2::Construct_curve_2 */

/*!


Functor to construct a `X_monotone_curve_2`. To enable caching the class 
is not default constructible and must be obtained via the function 
`construct_x_monotone_curve_2_object()`, which is a member of the traits. 

\cgalModels `Assignable`
\cgalModels `CopyConstructible`
\cgalModels `AdaptableBinaryFunction`
\cgalModels `AdaptableUnaryFunction`

*/
class Construct_x_monotone_curve_2 {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef AlgebraicKernel_d_1::Polynomial_1 Polynomial_1; 

/*!

*/ 
typedef AlgebraicKernel_d_1::Algebraic_real_1 Algebraic_real_1; 

/*!

*/ 
typedef Arr_rational_function_traits_2<AlgebraicKernel_d_1>::X_monotone_curve_2 result_type; 

/*!

*/ 
typedef Polynomial_1 argument_type; 

/*!

*/ 
typedef Polynomial_1 first_argument_type; 

/*!

*/ 
typedef Polynomial_1 second_argument_type; 

/// @} 

/// \name Operations 
/// @{

/*!
Constructs an \f$ x\f$-monotone curve supported by the polynomial function 
\f$ y = P(x)\f$. 
*/ 
X_monotone_curve_2 operator()(Polynomial_1 P) const; 

/*!
Constructs an \f$ x\f$-monotone curve supported by the polynomial function 
\f$ y = P(x)\f$. The function is defined over the interval \f$ [x,+\infty)\f$ if 
\f$ right\f$ is true and \f$ (-\infty,x]\f$ otherwise. 
*/ 
X_monotone_curve_2 operator()(Polynomial_1 P, 
const Algebraic_real_1& x, 
bool right) const; 

/*!
Constructs an \f$ x\f$-monotone curve supported by the polynomial function 
\f$ y = P(x)\f$. The function is defined over the interval \f$ [lower,upper]\f$. 
*/ 
X_monotone_curve_2 operator()(Polynomial_1 P, 
const Algebraic_real_1& lower, 
const Algebraic_real_1& upper); const 

/*!
Constructs an \f$ x\f$-monotone curve supported by the rational function 
\f$ y = P(x)/Q(x)\f$. 
\pre \f$ Q\f$ has no real roots. 
*/ 
X_monotone_curve_2 operator()(Polynomial_1 P, Polynomial_1 Q); const 

/*!
Constructs an \f$ x\f$-monotone curve supported by the rational function 
\f$ y = P(x)/Q(x)\f$. The function is defined over the interval \f$ I=[x,+\infty)\f$ 
if \f$ right\f$ is true and \f$ I=(-\infty,x]\f$ otherwise. 
\pre \f$ Q\f$ has no real roots in the interior of \f$ I\f$. 
*/ 
X_monotone_curve_2 operator()(Polynomial_1 P, Polynomial_1 Q, 
const Algebraic_real_1& x, 
bool right); const 

/*!
Constructs an \f$ x\f$-monotone curve supported by the rational function 
\f$ y = P(x)/Q(x)\f$. The function is defined over the interval \f$ I=[lower,upper]\f$. 
\pre \f$ Q\f$ has no real roots in the interior of \f$ I\f$. 
*/ 
X_monotone_curve_2 operator()(Polynomial_1 P, Polynomial_1 Q, 
const Algebraic_real_1& lower, 
const Algebraic_real_1& upper); const 

/*!
Constructs an \f$ x\f$-monotone curve supported by the polynomial function 
\f$ y = P(x)\f$, where the coefficients of \f$ P\f$ are given in the range 
`[begin,end)`. 
*/ 
template <typename InputIterator> 
X_monotone_curve_2 operator()(InputIterator begin, InputIterator end) const; 

/*!
Constructs an \f$ x\f$-monotone curve supported by the polynomial function 
\f$ y = P(x)\f$, where the coefficients of \f$ P\f$ are given in the range 
`[begin,end)`. The function is defined over the interval \f$ [x,+\infty)\f$ 
if \f$ right\f$ is true and \f$ (-\infty,x]\f$ otherwise. 
*/ 
template <typename InputIterator> 
X_monotone_curve_2 operator()(InputIterator begin, InputIterator end, 
const Algebraic_real_1& x, bool right) const; 

/*!
Constructs an \f$ x\f$-monotone curve supported by the polynomial function 
\f$ y = P(x)\f$, where the coefficients of \f$ P\f$ are given in the range 
`[begin,end)`. The function is defined over the interval 
\f$ [lower,upper]\f$. 
*/ 
template <typename InputIterator> 
X_monotone_curve_2 operator()(InputIterator begin, InputIterator end 
const Algebraic_real_1& lower, 
const Algebraic_real_1& upper); const 

/*!
Constructs an \f$ x\f$-monotone curve supported by the rational function 
\f$ y = P(x)/Q(x)\f$, where the coefficients of \f$ P\f$ and \f$ Q\f$ are given in the 
ranges `[begin_numer,end_numer)` and `[begin_denom,end_denom)`, 
respectively. 
\pre \f$ Q\f$ has no real roots. 
*/ 
template <typename InputIterator> 
X_monotone_curve_2 operator()(InputIterator begin_numer, InputIterator end_numer, 
InputIterator begin_denom, InputIterator end_denom); const 

/*!
Constructs an \f$ x\f$-monotone curve supported by the rational function 
\f$ y = P(x)/Q(x)\f$, where the coefficients of \f$ P\f$ and \f$ Q\f$ are given in the 
ranges `[begin_numer,end_numer)` and `[begin_denom,end_denom)`, 
respectively. The function is defined over the interval \f$ I=[x,+\infty)\f$ 
if \f$ right\f$ is true and \f$ I=(-\infty,x]\f$ otherwise. 
\pre \f$ Q\f$ has no real roots in the interior of \f$ I\f$. 
*/ 
template <typename InputIterator> 
X_monotone_curve_2 operator()(InputIterator begin_numer, InputIterator end_numer, 
InputIterator begin_denom, InputIterator end_denom, 
const Algebraic_real_1& x, bool right); const 

/*!
Constructs an \f$ x\f$-monotone curve supported by the rational function 
\f$ y = P(x)/Q(x)\f$, where the coefficients of \f$ P\f$ and \f$ Q\f$ are given in the 
ranges `[begin_numer,end_numer)` and `[begin_denom,end_denom)`, 
respectively. The function is defined over the interval \f$ I=[lower,upper]\f$. 
\pre \f$ Q\f$ has no real roots in the interior of \f$ I\f$. 
*/ 
template <typename InputIterator> 
X_monotone_curve_2 operator()(InputIterator begin_numer, InputIterator end_numer, 
InputIterator begin_denom, InputIterator end_denom, 
const Algebraic_real_1& lower, const Algebraic_real_1& upper); const 

/// @}

}; /* end Arr_rational_function_traits_2::Construct_x_monotone_curve_2 */

/*!


The `Curve_2` class nested within the traits is used 
to represent rational functions which may be restricted to a certain x-range. 

\cgalModels `ArrTraits::Curve_2`

*/
class Curve_2 {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef AlgebraicKernel_d_1::Polynomial_1 Polynomial_1; 

/*!

*/ 
typedef AlgebraicKernel_d_1::Algebraic_real_1 Algebraic_real_1; 

/// @} 

/// \name Operations 
/// @{

/*!
returns the numerator of the supporting rational function. 
*/ 
const Polynomial_1& numerator () const; 

/*!
returns the denominator of the supporting rational function. 
*/ 
const Polynomial_1& denominator () const; 

/*!
returns whether &ccedil;urve is continuous, namely whether it does not 
contains any vertical asymptotes in its interior. 
*/ 
bool is_continuous() const; 

/*!
returns whether the \f$ x\f$-coordinate of &ccedil;urve's left end is finite or 
whether it is \f$ \pm\infty\f$. 
*/ 
Arr_parameter_space left_parameter_space_in_x () const; 

/*!
returns whether the \f$ x\f$-coordinate of &ccedil;urve's right end is finite or 
whether it is \f$ \pm\infty\f$. 
*/ 
Arr_parameter_space right_parameter_space_in_x () const; 

/*!
returns the \f$ x\f$-coordinate of the left end. 
\pre left_boundary_in_x()==ARR_INTERIOR 
*/ 
Algebraic_real_1 left_x() const; 

/*!
returns the \f$ x\f$-coordinate of the right end. 
\pre right_boundary_in_x()==ARR_INTERIOR 
*/ 
Algebraic_real_1 right_x() const; 

/// @}

}; /* end Arr_rational_function_traits_2::Curve_2 */

/*!


\cgalModels `ArrTraits::Point_2`

*/
class Point_2 {
public:

/// \name Types 
/// @{

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

/// \name Operations 
/// @{

/*!
returns the numerator of the supporting rational function. 
*/ 
Polynomial_1 numerator () const; 

/*!
returns the denominator of the supporting rational function. 
*/ 
Polynomial_1 denominator () const; 

/*!
returns double-approximations of the x- and y-coordinates. 
*/ 
std::pair<double,double> to_double() const; 

/*!
returns the \f$ x\f$-coordinate of the point. 
*/ 
Algebraic_real_1 x() const; 

/*!
obtains the y-coordinates of the point. <B>Attention:</B> As described above, 
points are not stored by their y-coordinate in `Algebraic_real_1` 
representation. In fact, this representation must be computed on demand, and 
might become quite costly for points defined by high-degree polynomials. 
Therefore, it is recommended to avoid calls to this function as much as 
possible. 
*/ 
Algebraic_real_1 y() const; 

/*!
Computes a pair \f$ p\f$ approximating the \f$ x\f$-coordinate with 
respect to the given absolute precision \f$ a\f$. 
\post \f$ p.first \leq x \leq p.second \f$ 
\post \f$ p.second - p.first \leq2^{-a} \f$ 
*/ 
std::pair<Bound,Bound> approximate_absolute_x(int a) const; 

/*!
Computes a pair \f$ p\f$ approximating the \f$ y\f$-coordinate with 
respect to the given absolute precision \f$ a\f$. 
\post \f$ p.first \leq y \leq p.second \f$ 
\post \f$ p.second - p.first \leq2^{-a} \f$ 
*/ 
std::pair<Bound,Bound> approximate_absolute_y(int a) const; 

/*!
Computes a pair \f$ p\f$ approximating the \f$ x\f$-coordinate with 
respect to the given relative precision \f$ r\f$. 
\post \f$ p.first \leq x \leq p.second \f$ 
\post \f$ p.second - p.first \leq2^{-r}|x| \f$ 
*/ 
std::pair<Bound,Bound> approximate_relative_x(int r) const; 

/*!
Computes a pair \f$ p\f$ approximating the \f$ y\f$-coordinate with 
respect to the given relative precision \f$ r\f$. 
\post \f$ p.first \leq y \leq p.second \f$ 
\post \f$ p.second - p.first \leq2^{-r}|y| \f$ 
*/ 
std::pair<Bound,Bound> approximate_relative_y(int r) const; 

/// @}

}; /* end Arr_rational_function_traits_2::Point_2 */

/*!


The `X_monotone_curve_2` class nested within the traits is used 
to represent \f$ x\f$-monotone parts of rational functions. In particular, such an \f$ x\f$-monotone curve 
may not contain a vertical asymptote in its interior \f$ x\f$-range. 

\cgalModels `ArrTraits::XMonotoneCurve_2`

*/
class X_monotone_curve_2 {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef AlgebraicKernel_d_1::Polynomial_1 Polynomial_1; 

/*!

*/ 
typedef AlgebraicKernel_d_1::Algebraic_real_1 Algebraic_real_1; 

/*!

*/ 
typedef Arr_rational_function_traits_2<AlgebraicKernel_d_1>::Point_2 Point_2; 

/// @} 

/// \name Operations 
/// @{

/*!
returns the numerator of the supporting rational function. 
*/ 
const Polynomial_1& numerator () const; 

/*!
returns the denominator of the supporting rational function. 
*/ 
const Polynomial_1& denominator () const; 

/*!
returns whether the \f$ x\f$-coordinate of the source is finite or 
whether it is \f$ \pm\infty\f$. 
*/ 
Arr_parameter_space source_parameter_space_in_x () const; 

/*!
returns whether the \f$ y\f$-coordinate of the source is finite or 
whether it is \f$ \pm\infty\f$. 
*/ 
Arr_parameter_space source_parameter_space_in_y () const; 

/*!
returns the source point of the arc. 
\pre Both the \f$ x\f$- and \f$ y\f$-coordinates of the source point is finite. 
*/ 
const Point_2& source() const; 

/*!
returns the \f$ x\f$-coordinate of the source point. 
\pre The \f$ x\f$-coordinate of the source point is finite. 
*/ 
Algebraic_real_1 source_x() const; 

/*!
returns whether the \f$ x\f$-coordinate of the target is finite or 
whether it is \f$ \pm\infty\f$. 
*/ 
Arr_parameter_space target_parameter_space_in_x () const; 

/*!
returns whether the \f$ y\f$-coordinate of the target is finite or 
whether it is \f$ \pm\infty\f$. 
*/ 
Arr_parameter_space target_parameter_space_in_y () const; 

/*!
returns the target point of the arc. 
\pre Both the \f$ x\f$- and \f$ y\f$-coordinates of the target point is finite. 
*/ 
const Point_2& target() const; 

/*!
returns the \f$ x\f$-coordinate of the target point. 
\pre The \f$ x\f$-coordinate of the target point is finite. 
*/ 
Algebraic_real_1 target_x() const; 

/*!
returns whether the \f$ x\f$-coordinate of the left curve end is finite or 
whether it is \f$ \pm\infty\f$. 
*/ 
Arr_parameter_space left_parameter_space_in_x () const; 

/*!
returns whether the \f$ y\f$-coordinate of the left curve end is finite or 
whether it is \f$ \pm\infty\f$. 
*/ 
Arr_parameter_space left_parameter_space_in_y () const; 

/*!
returns the left point of the arc. 
\pre Both the \f$ x\f$- and \f$ y\f$-coordinates of the left point is finite. 
*/ 
const Point_2& left() const; 

/*!
returns the \f$ x\f$-coordinate of the left point. 
\pre The \f$ x\f$-coordinate of the left point is finite. 
*/ 
Algebraic_real_1 left_x() const; 

/*!
returns whether the \f$ x\f$-coordinate of the right curve end is finite or 
whether it is \f$ \pm\infty\f$. 
*/ 
Arr_parameter_space right_parameter_space_in_x () const; 

/*!
returns whether the \f$ y\f$-coordinate of the right curve end is finite or 
whether it is \f$ \pm\infty\f$. 
*/ 
Arr_parameter_space right_parameter_space_in_y () const; 

/*!
returns the right point of the arc. 
\pre Both the \f$ x\f$- and \f$ y\f$-coordinates of The right point is finite. 
*/ 
const Point_2& right() const; 

/*!
returns the \f$ x\f$-coordinate of the right point. 
\pre The \f$ x\f$-coordinate of the right point is finite. 
*/ 
Algebraic_real_1 right_x() const; 

/*!
returns whether the curve is oriented from left to right. 
*/ 
bool is_left_to_right () const; 

/// @}

}; /* end Arr_rational_function_traits_2::X_monotone_curve_2 */


}; /* end Arr_rational_function_traits_2 */
} /* end namespace CGAL */
