/*!
\ingroup PkgConeBasedSpannersConcepts
\cgalConcept

The concept `ComputeConeBoundaries_2` describes the set of requirements for the functor that 
computes the directions of cone boundaries with a given cone number and a given initial direction 
either exactly or inexactly.

\cgalHasModel `CGAL::Compute_cone_boundaries_2`

*/
template <typename Kernel_>
class ComputeConeBoundaries_2 {
public:

/// \name Types 
/// @{

/*! The CGAL kernel type used by the functor. If this parameter is
	`CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt`,
    the cone boundaries will be computed exactly; otherwise, inexactly using an approximate PI=3.1415...
	This kernel type also decides other types such as Kernel_::Point_2, Kernel_::Direction_2, etc.
*/
typedef Kernel_                          kernel_type;

/// @}


/// \name Operator 
/// @{

/*! \brief The operator().
 *
 *  Compute the directions of cone boundaries with a given
 *  cone number and a given initial direction. The results are returned by the reference
 *  argument: vector \p rays.
 *
 *  \param[in] cone_number The number of cones
 *  \param[in] initial_direction The direction of the first ray
 *  \param[out] rays  Storing the results, a vector of directions. It should contain no
 *                    elements when passed to this operator.
 */
void operator()(const unsigned int cone_number,
				Direction_2& initial_direction,
				std::vector<Direction_2>& rays);

/// @} 

}; /* end ComputeConeBoundaries_2 */

