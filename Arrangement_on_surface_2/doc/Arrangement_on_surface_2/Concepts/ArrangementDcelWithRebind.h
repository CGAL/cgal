
/*!
\ingroup PkgArrangementOnSurface2ConceptsDCEL
\cgalConcept

The concept `ArrangementDcelWithRebind` refines the `ArrangementDcel` concept by adding
a policy clone idiom in form of a rebind struct-template.

Instantiate a dcel class with many different possible types without ad-hoc limitations on type of the dcel classes.

\cgalRefines{ArrangementDcel}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Arr_default_dcel<Traits>}
\cgalHasModels{CGAL::Arr_face_extended_dcel<Traits,FData,V,H,F>}
\cgalHasModels{CGAL::Arr_extended_dcel<Traits,VData,HData,FData,V,H,F>}
\cgalHasModelsEnd

*/

class ArrangementDcelWithRebind {
public:

/// \name Types
/// @{

/*!
allows the instantiation of a model of the base concept
`ArrangementDcel` with a different possible geometry-traits
class without ad-hoc limitations on it.

Following the standard clone policy, the rebind struct-template must
have a nested type named `other` that defines the type of the
model replica.
*/
typedef unspecified_type template <class T> rebind;

/// @}

/// \name Creation
/// @{

/*!
constructs an empty \dcel with one unbounded face.
*/
Arr_dcel();

/// @}

}; /* end ArrangementDcelWithRebind */

