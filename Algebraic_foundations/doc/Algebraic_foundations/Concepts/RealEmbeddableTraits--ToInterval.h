
namespace RealEmbeddableTraits_ {

/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

`AdaptableUnaryFunction` computes for a given real embeddable
number \f$ x\f$ a double interval containing \f$ x\f$.
This interval is represented by `std::pair<double,double>`.

\cgalRefines{AdaptableUnaryFunction}

\sa `RealEmbeddableTraits`

*/

class ToInterval {
public:

/// \name Types
/// @{

/*!
The result type.
*/
typedef std::pair<double,double> result_type;

/*!
Is `RealEmbeddableTraits::Type`.
*/
typedef unspecified_type argument_type;

/// @}

/// \name Operations
/// @{

/*!
computes a double interval containing \f$ x\f$.
*/
result_type operator()(argument_type x);

/// @}

}; /* end ToInterval */

} /* end of namespace RealEmbeddableTraits_ */
