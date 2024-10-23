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

/*!  must be `true` if  `Kernel` has a nested type `Has_filtered_predicates_tag`  with `value == true`
*/
    static constexpr bool is_filtered = unspecified;

    /*!  Must be `true` if the number type of `Kernel` is a floating point number type*/
    static constexpr bool  is_floating_point = unspecified;
/// \name Types
/// @{

/*! The kernel type. If this type has a nested type `Has_filtered_predicates_tag`  with `value == true`,
*/
using Kernel = unspecified;

/*! The point type of `Kernel` corresponding to `dimension`
*/
using Point = unspecified_type;

/*! The number type of the approximate kernel.
*/
using distance_t = Approximate_kernel::FT;

/*! The approximate kernel
*/
using Approximate_kernel = unspecified_type;


/*! The point type of the approximate kernel corresponding to `dimension`.
*/
using Approximate_point = unspecified_type;

/*! A functor of the approximate kernel for filtered points
*/
using Construct_bbox = unspecified_type;

/*! A functor of the approximate kernel for two filtered points
*/
using Squared_distance = unspecified_type;


/*! A functor of the approximate kernel for two filtered points
*/
using Difference_of_points = unspecified_type;


/*! A functor of the approximate kernel the return type of `Difference_of_points`
*/
using Scaled_vector = unspecified_type;


/*! A functor of the approximate kernel the return type of `Scaled_vector`
*/
using Translated_point = unspecified_type;

/*! The exact kernel
*/
using Exact_kernel = unspecified_type;

/*! The point type of the exact kernel corresponding to `dimension`.
The  point type must have `operator[]` returning a number type which can be used as first template parameter of `Sqrt_extension`.
*/
using Exact_point = unspecified_type;

/*! A converter for points from `Kernel` to  `Approximate_kernel`
*/
using K2A = unspecified_type;


/*! A converter for points from `Approximate_kernel` to  `Exact_kernel`
*/
using A2E = unspecified_type;
/// @}
};
