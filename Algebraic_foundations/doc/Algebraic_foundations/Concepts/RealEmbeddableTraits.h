
/*!
\ingroup PkgAlgebraicFoundationsRealEmbeddableConcepts
\cgalConcept

A model of `RealEmbeddableTraits` is associated to a number type 
`Type` and reflects the properties of this type with respect 
to the concept `RealEmbeddable`. 

\cgalHasModel `CGAL::Real_embeddable_traits<T>`
*/

class RealEmbeddableTraits {
public:

/// \name Types 
/// A model of `RealEmbeddableTraits` is supposed to provide:
/// @{

/*! 
The associated number type. 
*/ 
typedef Hidden_type Type; 

/*! 
Tag indicating whether the associated type is real embeddable. 

This is either \link Tag_true `Tag_true`\endlink or \link Tag_false `Tag_false`\endlink. 
*/ 
typedef Hidden_type Is_real_embeddable; 

/*! 
This type specifies the return type of the predicates provided 
by this traits. The type must be convertible to `bool` and 
typically the type indeed maps to `bool`. However, there are also 
cases such as interval arithmetic, in which it is `Uncertain<bool>` 
or some similar type. 

*/ 
typedef Hidden_type Boolean; 

/*! 
This type specifies the return type of the `Sgn` functor. 
The type must be convertible to `CGAL::Sign` and 
typically the type indeed maps to `CGAL::Sign`. However, there are also 
cases such as interval arithmetic, in which it is `Uncertain<CGAL::Sign>` 
or some similar type. 

*/ 
typedef Hidden_type Sign; 

/*! 
This type specifies the return type of the `Compare` functor. 
The type must be convertible to `CGAL::Comparison_result` and 
typically the type indeed maps to `CGAL::Comparison_result`. However, there are also 
cases such as interval arithmetic, in which it is `Uncertain<CGAL::Comparison_result>` 
or some similar type. 

*/ 
typedef Hidden_type Comparison_result; 

/// @} 

/// \name Functors 
/// In case the associated type is `RealEmbeddable` all functors are
/// provided. In case a functor is not provided, it is set to
/// `CGAL::Null_functor`.
/// @{

/*! 

A model of `RealEmbeddableTraits::IsZero` 
In case `Type` is also model of `IntegralDomainWithoutDivision` 
this is a model of `AlgebraicStructureTraits::IsZero`. 
*/ 
typedef Hidden_type Is_zero; 

/*! 
A model of `RealEmbeddableTraits::Abs` 
*/ 
typedef Hidden_type Abs; 

/*! 
A model of `RealEmbeddableTraits::Sgn` 
*/ 
typedef Hidden_type Sgn; 

/*! 
A model of `RealEmbeddableTraits::IsPositive` 
*/ 
typedef Hidden_type Is_positive; 

/*! 
A model of `RealEmbeddableTraits::IsNegative` 
*/ 
typedef Hidden_type Is_negative; 

/*! 
A model of `RealEmbeddableTraits::Compare` 
*/ 
typedef Hidden_type Compare; 

/*! 
A model of `RealEmbeddableTraits::ToDouble` 
*/ 
typedef Hidden_type To_double; 

/*! 
A model of `RealEmbeddableTraits::ToInterval` 
*/ 
typedef Hidden_type To_interval; 

/// @}

}; /* end RealEmbeddableTraits */

