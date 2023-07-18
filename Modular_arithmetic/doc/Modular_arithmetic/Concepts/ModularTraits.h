
/*!
\ingroup PkgModularArithmeticConcepts
\cgalConcept

A model of `ModularTraits` is associated to a specific `Type`.
In case this associated type is a model of `Modularizable`, this is indicated by the
Boolean tag `ModularTraits::Is_modularizable`. The mapping into the `Residue_type` is
provided by the functor `ModularTraits::Modular_image`.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Modular_traits<T>}
\cgalHasModelsEnd

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
typedef unspecified_type Type ;

/*!

Tag indicating whether the associated type is modularizable.

This is either `CGAL::Tag_true` or `CGAL::Tag_false`.
*/
typedef unspecified_type Is_modularizable;

/*!

The type of the modular image.

In case the associated type is not a model of `Modularizable` this is undefined.

*/
typedef unspecified_type Residue_type;

/// @}

/// \name Functors
/// In case the associated type is a model of `Modularizable` all functors are
/// provided. In case a functor is not provided, it is set to
/// `CGAL::Null_functor`.
/// @{

/*!
A model of `ModularTraits::ModularImage`
*/
typedef unspecified_type Modular_image;

/*!
A model of `ModularTraits::ModularImageRepresentative`
*/
typedef unspecified_type Modular_image_representative;

/// @}

}; /* end ModularTraits */

