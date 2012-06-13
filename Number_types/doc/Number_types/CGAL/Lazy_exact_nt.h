namespace CGAL {
/// \addtogroup PkgNumberTypes Number Types
/// @{

 
///  An object of the class `Lazy_exact_nt<NT>` is able to represent any 
///  real embeddable number which `NT` is able to represent.
///  The idea is that `Lazy_exact_nt<NT>` works exactly like `NT`, except
///  that it is expected to be faster because it tries to only compute an 
///  approximation of the value, and only refers to `NT` when needed.  
///  The goal is to speed up exact computations done by any exact but slow 
///  number type `NT`.
///
///  `NT` must be a model of concept `RealEmbeddable`. 
///  `NT` must be at least model of concept `IntegralDomainWithoutDivision`.
///
///  Note that some filtering mechanism is available at the predicate level
///  using `Filtered_predicate` and `Filtered_kernel`.
///  
///  \models ::IntegralDomainWithoutDivision same as  `NT`
///  \models ::RealEmbeddable
///  \models ::Fraction if `NT` is a ::Fraction
///
///  Example 
///  -------------- 
///  
///  \code{.cpp}
///  #include <CGAL/Cartesian.h>
///  #include <CGAL/MP_Float.h>
///  #include <CGAL/Lazy_exact_nt.h>
///  #include <CGAL/Quotient.h>
///  
///  typedef CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> > NT;
///  typedef CGAL::Cartesian<NT> K;
///  \endcode
///   
template< class NT >
class Lazy_exact_nt {
public:

/// \name Creation
/// @{
/*!
 introduces an uninitialized variable `m`.
*/
Lazy_exact_nt();
/// @}

/// \name Creation
/// @{
/*!
 introduces the integral value `i`.
*/
  Lazy_exact_nt(int i);
/// @}

/// \name Creation
/// @{
/*!
 introduces the floating point value `d` (works only if `NT` has a constructor from a double too).
*/
Lazy_exact_nt(double d);
/// @}

/// \name Creation
/// @{
/*!
 introduces the value `n`.
*/
  Lazy_exact_nt(NT n);
/// @}

/// \name Creation
/// @{
/*!
 introduces the value `n`. `NT1` needs to be convertible to `NT` (and this conversion will only be done if necessary).
*/
  template <class NT1> Lazy_exact_nt(Lazy_exact_nt<NT1> n);
/// @}

/// \name Operations
/// @{
/*!
 returns the corresponding NT value.
*/
NT exact();
/// @}

/// \name Operations
/// @{
/*!
 returns an interval containing the exact value.
*/
Interval_nt<true> approx();
/// @}

/// \name Operations
/// @{
/*!
 returns an interval containing the  exact value.
*/
Interval_nt<false> interval();
/// @}

/// \name Operations
/// @{
/*!
 specifies the relative precision that `to_double()` has to fulfill. The default value is \f$10^-5\f$.  \pre d>0 and d<1.
*/
static void set_relative_precision_of_to_double(double d);
/// @}

/// \name Operations
/// @{
/*!
 returns the relative precision that `to_double()` currently fulfills.
*/
static double get_relative_precision_of_to_double();
/// @}

}; /* class Lazy_exact_nt */

/// \relates Lazy_exact_nt
/// writes `m` to ostream `out` in an interval format.
std::ostream& operator<<(std::ostream& out,
                                     const Lazy_exact_nt<NT>& m);

/// \relates Lazy_exact_nt
/// reads a `NT` from `in`, then converts it to a `Lazy_exact_nt<NT>`.
std::istream& operator>>(std::istream& in, Lazy_exact_nt<NT>& m);

}; /* class Lazy_exact_nt */

/// @}
} // namespace CGAL

                   
  

