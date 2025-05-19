
/*!
\ingroup PkgModularArithmeticConcepts
\cgalConcept

An algebraic structure is called `Modularizable`, if there is a suitable mapping
into an algebraic structure which is based on the type `CGAL::Residue`.
For scalar types, e.g. Integers, this mapping is just the canonical homomorphism
into the type `CGAL::Residue` with respect to the current prime.
For compound types, e.g. Polynomials,
the mapping is applied to the coefficients of the compound type.

The mapping is provided via `CGAL::Modular_traits<Modularizable>`,
being a model of `ModularTraits`.

Note that types representing rationals, or types which do have some notion
of denominator, are not `Modularizable`.
This is due to the fact that the denominator may be zero modulo the prime,
which can not be represented.

\cgalHasModelsBegin
\cgalHasModels{int}
\cgalHasModels{long}
\cgalHasModels{CORE::BigInt}
\cgalHasModels{CGAL::Gmpz}
\cgalHasModels{leda_integer}
\cgalHasModels{mpz_class}
\cgalHasModelsEnd

The following types are `Modularizable` iff their template arguments are.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Lazy_exact_nt<NT>}
\cgalHasModels{CGAL::Sqrt_extension<NT,ROOT>}
\cgalHasModels{CGAL::Polynomial<Coeff>}
\cgalHasModelsEnd

\sa `CGAL::Residue`
\sa `CGAL::Modular_traits<T>`

*/

class Modularizable {
public:

}; /* end Modularizable */

