/*!
\ingroup PkgFrechetDistanceConcepts
\cgalConcept

The concept `FrechetDistanceTraits` defines the requirements ...


\cgalHasModelsBegin
\cgalHasModels{CGAL::Frechet_distance_traits_2}
\cgalHasModels{CGAL::Frechet_distance_traits_3}
\cgalHasModels{CGAL::Frechet_distance_traits_d}
\cgalHasModelsEnd
*/

class FrechetDistanceTraits {

    public:
/*!  2 or 3 */
    const int dimension;

/// \name Types
/// @{

/*! The kernel type. If this type has a nested type `Has_filtered_predicates_tag`  with `value == true`,
    it must have a nested type `Exact_kernel`, and a nested type `C2E` with an  `operator()` that converts
    a point of `Kernel` to a point of `Exact_kernel`. Otherwise, it must have a nested type `FT` for
    which an overload of `to_double()` exists.
*/
using Kernel = unspecified;

/*! The point type of `Kernel` corresponding to `dimension`
*/
using Point = unspecified_type;

/// @}
};