
/*!
\ingroup PkgModularArithmeticConcepts
\cgalconcept

A model of `ModularTraits` is associated to a specific `Type`. 
In case this associated type is `Modularizable`, this is indicated by the 
boolean tag `Is_modularizable`. The mapping into the `Residue_type` is 
provided by the functor `Modular_image`. 

\hasModel CGAL::Modular_traits<T> 

\sa `CGAL::Residue` 
\sa `Modularizable` 

*/

class ModularTraits {
public:

/// \name Types 
/// A model of `ModularTraits` is supposed to provide:
/// @{

/*! 
The associated type. 
*/ 
typedef Hidden_type Type ; 

/*! 

Tag indicating whether the associated type is modularizable. 

This is either `CGAL::Tag_true` or `CGAL::Tag_false`. 
*/ 
typedef Hidden_type Is_modularizable; 

/*! 

The type of the modular image. 

In case the type is not `Modularizable` this is undefined. 

*/ 
typedef Hidden_type Residue_type; 

/// @} 

/// \name Functors 
/// In case the associated type is `Modularizable` all functors are
/// provided. In case a functor is not provided, it is set to
/// `CGAL::Null_functor`.
/// @{

/*! 
A model of `ModularTraits::ModularImage` 
*/ 
typedef Hidden_type Modular_image; 

/*! 
A model of `ModularTraits::ModularImageRepresentative` 
*/ 
typedef Hidden_type Modular_image_representative; 

/// @}

}; /* end ModularTraits */

