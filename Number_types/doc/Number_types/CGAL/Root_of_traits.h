///
namespace CGAL {
/// \addtogroup PkgNumberTypes Number Types
/// @{
///
 
///  
///  The function `compute_roots_of_2` solves a univariate polynomial as it is defined by the 
///  coefficients given to the function. The solutions are written into the given
///  `OutputIterator`.
///  
///  \sa ::RootOf_2
///  \sa `CGAL::Root_of_traits<RT>`
///  \sa `CGAL::make_root_of_2<RT>`
///  \sa `CGAL::make_sqrt<RT>`
///  \sa `CGAL::Sqrt_extension<NT,ROOT>`
/// Writes the real roots of the polynomial \f$aX^2+bX+c\f$ into \f$oit\f$ in ascending order.  `OutputIterator` is required to accept `Root_of_traits::Root_of_2`.  Multiplicities are not reported.  \pre `RT` is an `IntegralDomainWithoutDivision`. \pre \f$a \neq 0\f$ or \f$b \neq 0\f$. 
template <typename RT, typename OutputIterator>
            OutputIterator 
            compute_roots_of_2(const RT& a, const RT& b, const RT& c, OutputIterator oit);
/// @}
} 

///
namespace CGAL {
/// \addtogroup PkgNumberTypes Number Types
/// @{
///
 
///  
///  The function `make_root_of_2` constructs an algebraic number of degree 2 over a 
///  ring number type.
///  
///  \sa `RootOf_2`
///  \sa `CGAL::Root_of_traits<RT>`
///  \sa `CGAL::make_sqrt<RT>`
///  \sa `CGAL::compute_roots_of_2<RT,OutputIterator>`
///  \sa `CGAL::Sqrt_extension<NT,ROOT>`
/// Returns the smallest real root of the polynomial \f$aX^2+bX+c\f$ if \f$s\f$ is true, and the largest root is \f$s\f$ is false. \pre `RT` is an `IntegralDomainWithoutDivision`. \pre The polynomial has at least one real root.
template <typename RT>
            Root_of_traits<RT>::Root_of_2
            make_root_of_2(const RT& a, const RT& b, const RT& c, bool s);
/// @}
} 

///
namespace CGAL {
/// \addtogroup PkgNumberTypes Number Types
/// @{
///
 
///  
///  The function `make_root_of_2` constructs an algebraic number of degree 2 over a 
///  ring number type.
///  
///  \sa `RootOf_2`
///  \sa `CGAL::Root_of_traits<RT>`
///  \sa `CGAL::make_sqrt<RT>`
///  \sa `CGAL::compute_roots_of_2<RT,OutputIterator>`
///  \sa `CGAL::Sqrt_extension<NT,ROOT>`
/// Constructs the number \f$&alpha;+ &beta; \sqrt{&gamma;} \f$.  \pre `RT` is an `IntegralDomainWithoutDivision`. \pre \f$&gamma;  \geq 0\f$
template <typename RT>
Root_of_traits<RT>::Root_of_2
make_root_of_2(RT alpha, RT beta, RT gamma);
/// @}
} 

///
namespace CGAL {
/// \addtogroup PkgNumberTypes Number Types
/// @{
///
 
///  
///  The function `make_sqrt` constructs a square root of a given value of type \f$RT\f$. 
///  Depending on the type \f$RT\f$ the square root may be returned in a new type that 
///  can represent algebraic extensions of degree \f$2\f$.
///  
///  \sa `RootOf_2`
///  \sa `CGAL::make_root_of_2<RT>`
///  \sa `CGAL::Root_of_traits<RT>`
/// Returns \f$ \sqrt{x} .\f$ \pre `RT` is a `RealEmbeddable` `IntegralDomain`. \pre  \f$x  \leq 0 \f$ 
template <typename RT> Root_of_traits<RT>::Root_of_2 make_sqrt(const RT& x);
/// @}
} 

namespace CGAL {
/// \addtogroup PkgNumberTypes Number Types
/// @{

 
///  
///  For a `RealEmbeddable` `IntegralDomain` `RT`, the class template 
///  `Root_of_traits<RT>` associates a type `Root_of_2`, which represents 
///  algebraic numbers of degree 2 over `RT`. Moreover, the class provides 
///  `Root_of_1`, which represents the quotient field of `RT`.
///  
///  
///  
///  \sa `RootOf_2`
///  \sa `CGAL::compute_roots_of_2<RT,OutputIterator>`
///  \sa `CGAL::make_root_of_2<RT>`
template< class RT >
class Root_of_traits {
public:

/// \name Types
/// @{
/*!
 A `RealEmbeddable` `Field` representing the quotient field of `RT`.
*/
typedef Hidden_type Root_of_1;
/// @}

/// \name Types
/// @{
/*!
 Model of `RootOf_2`.
*/
typedef Hidden_type Root_of_2;
/// @}

}; /* class Root_of_traits */
/// @}
} // namespace CGAL

                   
  

