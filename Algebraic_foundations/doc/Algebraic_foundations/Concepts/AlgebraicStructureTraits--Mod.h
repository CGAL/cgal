
namespace AlgebraicStructureTraits_ {

/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

`AdaptableBinaryFunction` computes the remainder of division with remainder.

\cgalRefines{AdaptableBinaryFunction}

\sa `AlgebraicStructureTraits`
\sa `AlgebraicStructureTraits_::Div`
\sa `AlgebraicStructureTraits_::DivMod`

*/

class Mod {
public:

/// \name Types
/// @{

/*!
Is `AlgebraicStructureTraits::Type`.
*/
typedef unspecified_type result_type;

/*!
Is `AlgebraicStructureTraits::Type`.
*/
typedef unspecified_type first_argument;

/*!
Is `AlgebraicStructureTraits::Type`.
*/
typedef unspecified_type second_argument;

/// @}

/// \name Operations
/// @{

/*!

*/
result_type operator()(first_argument_type x,
second_argument_type y);

/*!
This operator is defined if `NT1` and `NT2` are `ExplicitInteroperable`
with coercion type `AlgebraicStructureTraits::Type`.
*/
template <class NT1, class NT2> result_type operator()(NT1 x, NT2 y);

/// @}

}; /* end Mod */

} /* end of namespace AlgebraicStructureTraits_ */
