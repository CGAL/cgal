
namespace AlgebraicStructureTraits_{

/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

`AdaptableBinaryFunction` providing an integral division. 

Integral division (a.k.a. exact division or division without remainder) maps 
ring elements \f$ (x,y)\f$ to ring element \f$ z\f$ such that \f$ x = yz\f$ if such a \f$ z\f$ 
exists (i.e.\ if \f$ x\f$ is divisible by \f$ y\f$). Otherwise the effect of invoking 
this operation is undefined. Since the ring represented is an integral domain, 
\f$ z\f$ is uniquely defined if it exists. 

\cgalRefines `AdaptableBinaryFunction` 

\sa `AlgebraicStructureTraits`
\sa `AlgebraicStructureTraits_::Divides`

*/

class IntegralDivision {
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
returns \f$ x/y\f$, this is an integral division. 
*/ 
result_type operator()(first_argument_type x, 
second_argument_type y); 

/*!
This operator is defined if `NT1` and `NT2` are `ExplicitInteroperable` 
with coercion type `AlgebraicStructureTraits::Type`. 
*/ 
template <class NT1, class NT2> result_type operator()(NT1 x, NT2 y); 

/// @}

}; /* end IntegralDivision */

} /* end of namespace AlgebraicStructureTraits_ */


