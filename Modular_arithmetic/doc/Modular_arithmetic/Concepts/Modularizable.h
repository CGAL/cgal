/// \ingroup PkgModularArithmeticConcepts
 /// Definition 
/// --------------
/// 
/// An algebraic structure is called `Modularizable`, if there is a suitable mapping 
/// into an algebraic structure which is based on the type `CGAL::Residue`. 
/// For scalar types, e.g. Integers, this mapping is just the canonical homomorphism 
/// into the type `CGAL::Residue` with respect to the current prime. 
/// For compound types, e.g. Polynomials, 
/// the mapping is applied to the coefficients of the compound type.
/// The mapping is provided via `CGAL::Modular_traits<Modularizable>`, 
/// being a model of `ModularTraits`.
/// Note that types representing rationals, or types which do have some notion 
/// of denominator, are not `Modularizable`. 
/// This is due to the fact that the denominator may be zero modulo the prime, 
/// which can not be represented.
/// 
/// 
/// Has Models 
/// --------------
/// 
/// `int`
/// `long`
/// `CORE::BigInt`
/// `CGAL::Gmpz`
/// `leda::integer`
/// `mpz_class`
/// The following types are `Modularizable` iff their template arguments are. 
/// `CGAL::Lazy_exact_nt<NT>` 
/// `CGAL::Sqrt_extension<NT,ROOT>`
/// `CGAL::Polynomial<Coeff>`
/// See Also 
/// --------------
/// `CGAL::Residue`
/// `CGAL::Modular_traits<T>`
///  
/// }; /* concept Modularizable */
class Modularizable {
public:

}; /* concept Modularizable */
///

                   
  

