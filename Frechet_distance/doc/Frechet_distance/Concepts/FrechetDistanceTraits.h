/*!
\ingroup PkgFrechetDistanceConcepts
\cgalConcept

The concept `FrechetDistanceTraits` defines the requirements of the
first template parameter of the functions `CGAL::is_Frechet_distance_larger()`
and `CGAL::approximate_Frechet_distance()`.


\cgalHasModelsBegin
\cgalHasModels{CGAL::Frechet_distance_traits_2}
\cgalHasModels{CGAL::Frechet_distance_traits_3}
\cgalHasModels{CGAL::Frechet_distance_traits_d}
\cgalHasModelsEnd
*/

class FrechetDistanceTraits {

    public:
/*!  a fixed dimension >= 2 */
    const int dimension;

/*!  a fixed dimension >= 2 */
    static constexpr bool is_filtered;
    /*!  a fixed dimension >= 2 */
    static constexpr bool  is_floating_point;
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

/*! The number type of the filtered kernel. If
*/
using distance_t = Filtered_kernel::FT;

/*! The filtered kernel
*/
using Filtered_kernel = unspecified_type;


/*! The point type of the filtered kernel corresponding to `dimension`.
*/
using Filtered_point = unspecified_type;

/*! A functor of the filtered kernel for filtered points
*/
using Construct_bbox = unspecified_type;

/*! A functor of the filtered kernel for two filtered points
*/
using Squared_distance = unspecified_type;


/*! A functor of the filtered kernel for two filtered points
*/
using Difference_of_points = unspecified_type;


/*! A functor of the filtered kernel the return type of `Difference_of_points`
*/
using Scaled_vector = unspecified_type;


/*! A functor of the filtered kernel the return type of `Scaled_vector`
*/
using Translated_point = unspecified_type;

/*! The exact kernel
*/
using Exact_kernel = unspecified_type;

/*! The point type of the exact kernel corresponding to `dimension`.
The  point type must have `operator[]` returning a number type which can be used as first template parameter of `Sqrt_extension`.
*/
using Exact_point = unspecified_type;

/*! A converter for points from `Kernel` to  `Filtered_kernel`
*/
using K2F = unspecified_type;


/*! A converter for points from `Filtered_kernel` to  `Exact_kernel`
*/
using F2E = unspecified_type;
/// @}
};
