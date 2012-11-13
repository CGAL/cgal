
/*!
\ingroup PkgAlgebraicFoundationsRealEmbeddableConcepts
\cgalConcept

A model of this concepts represents numbers that are embeddable on the real 
axis. The type obeys the algebraic structure and compares two values according 
to the total order of the real numbers. 

Moreover, `CGAL::Real_embeddable_traits< RealEmbeddable >` is a model of 
`RealEmbeddableTraits` 

with: 

- `CGAL::Real_embeddable_traits< RealEmbeddable >::Is_real_embeddable` set to `Tag_true` 

and functors : 

- `CGAL::Real_embeddable_traits< RealEmbeddable >::Is_zero` 

- `CGAL::Real_embeddable_traits< RealEmbeddable >::Abs` 

- `CGAL::Real_embeddable_traits< RealEmbeddable >::Sgn` 

- `CGAL::Real_embeddable_traits< RealEmbeddable >::Is_positive` 

- `CGAL::Real_embeddable_traits< RealEmbeddable >::Is_negative` 

- `CGAL::Real_embeddable_traits< RealEmbeddable >::Compare` 

- `CGAL::Real_embeddable_traits< RealEmbeddable >::To_double` 

- `CGAL::Real_embeddable_traits< RealEmbeddable >::To_interval` 

Remark: 

If a number type is a model of both `IntegralDomainWithoutDivision` and 
`RealEmbeddable`, it follows that the ring represented by such a number type 
is a sub-ring of the real numbers and hence has characteristic zero. 

\cgalRefines `Equality` Comparable 
\cgalRefines `LessThanComparable` 

\sa ::RealEmbeddableTraits 

*/

class RealEmbeddable {
public:

/// \name Operations 
/// @{

/*! 

*/ 
bool operator==(const RealEmbeddable &a, 
const RealEmbeddable &b); 


/*! 

*/ 
bool operator!=(const RealEmbeddable &a, 
const RealEmbeddable &b); 

/*! 

*/ 
bool operator< (const RealEmbeddable &a, 
const RealEmbeddable &b); 

/*! 

*/ 
bool operator<=(const RealEmbeddable &a, 
const RealEmbeddable &b); 

/*! 


*/ 
bool operator> (const RealEmbeddable &a, 
const RealEmbeddable &b); 

/*! 

\relates RealEmbeddable 
*/ 
bool operator>=(const RealEmbeddable &a, 
const RealEmbeddable &b); 

/// @}

}; /* end RealEmbeddable */

