namespace CGAL {
/// \addtogroup PkgNumberTypes Number Types
/// @{

 
///  An object of the class `Quotient` is an element of the 
///  field of quotients of the integral domain type `NT`.
///  If `NT` behaves like an integer, `Quotient`
///  behaves like a rational number. 
///
///  \leda's class `rational` (see Section  \ref ledant )
///  has been the basis for `Quotient`.
///  A `Quotient` `q` is represented as a pair of 
///  `NT`s, representing numerator and denominator.
///
///  - `NT` must be at least model of concept `IntegralDomainWithoutDivision`.
///  - `NT` must be a model of concept `RealEmbeddable`. 
///  
///  \models ::Field
///  \models ::RealEmbeddable
///  \models ::Fraction
template< class NT >
class Quotient {
public:
/// \name Creation
/// @{
/*!
 introduces an uninitialized variable `q`.
*/
Quotient();
/// @}

/// \name Creation
/// @{
/*!
 introduces the quotient `t/1`. NT needs to have a constructor from T.
*/
template <class T> Quotient(const T& t);
/// @}

/// \name Creation
/// @{
/*!
 introduces the quotient `NT(t.numerator())/NT(t.denominator())`. NT needs to have a constructor from T.
*/
template <class T> Quotient(const Quotient<T>& t);
/// @}

/// \name Creation
/// @{
/*!
 introduces the quotient `n/d`.             \pre \f$d  \neq 0\f$.         
*/
  Quotient(const NT& n, const NT& d);
/// @}


/// \name Operations
/// There are two access functions, namely to the
/// numerator and the denominator of a quotient.
/// Note that these values are not uniquely defined. 
/// It is guaranteed that `q.numerator()` and 
/// `q.denominator()` return values `nt_num` and
/// `nt_den` such that `q = nt_num/nt_den`, only
/// if  `q.numerator()` and `q.denominator()` are called
/// consecutively wrt `q`, i.e. `q` is not involved in 
/// any other operation between these calls.
///
/// @{

/*!
  returns a numerator of `q`.
*/
NT numerator() const;

/*!
  returns a denominator of `q`.
*/
NT denominator() const;

/// @}
 
}; /* class Quotient */

/// The stream operations are available as well. 
/// They assume that corresponding stream operators for type `NT` exist.
/// \relates Quotient
/// writes `q` to ostream `out` in format ``<TT>n/d</TT>'', where `n`\f$==\f$`q.numerator()` and `d`\f$==\f$`q.denominator()`.
std::ostream& operator<<(std::ostream& out, const Quotient<NT>& q);

///
/// The stream operations are available as well. 
/// They assume that corresponding stream operators for type `NT` exist.
/// \relates Quotient
/// reads `q` from istream `in`. Expected format is `n/d`, where `n` and `d` are of type `NT`. A single `n` which is not followed by a `/`  is also accepted and interpreted as `n/1`.
std::istream& operator>>(std::istream& in, Quotient<NT>& q);

/// \relates Quotient
/// returns some double approximation to `q`.
/// Provided to fulfill the \cgal requirements on number types.
double to_double(const Quotient<NT>& q);

/// \relates Quotient
/// returns true, if numerator and denominator are valid.
/// Provided to fulfill the \cgal requirements on number types.
bool  is_valid(const Quotient<NT>& q);

/// \relates Quotient
/// returns true, if numerator and denominator are finite.
/// Provided to fulfill the \cgal requirements on number types.
bool  is_finite(const Quotient<NT>& q);

/// \relates Quotient
/// returns the square root of `q`.  This is supported if and only if         `NT` supports the square root as well.
/// Provided to fulfill the \cgal requirements on number types.
Quotient<NT>  sqrt(const Quotient<NT>& q);

/// @}
} // namespace CGAL

                   
  

