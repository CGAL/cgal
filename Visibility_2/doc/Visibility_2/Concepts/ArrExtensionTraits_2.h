/*!
\ingroup PkgVisibility_2Concepts
\cgalConcept

All visibility polgyon algorithms provided in \cgal are parameterized with a traits class 'Traits', which defines the extension of Arrangement_2 the output will have.

\cgalHasModel `CGAL::Arr_extension_default_traits_2`

\sa `Visibility_2`

*/
class ArrExtensionTraits_2 {
public :

/// \name Creation
/// @{
/*!
copy creator
*/
ArrExtensionTraits_2 (const ArrExtensionTraits_2 & Traits);

/// @}

}

