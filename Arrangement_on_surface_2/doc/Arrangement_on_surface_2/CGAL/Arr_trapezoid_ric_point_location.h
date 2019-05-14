
namespace CGAL {

/*!
\ingroup PkgArrangement2PointLocation

\anchor arr_reftrap_pl 

The `Arr_trapezoid_ric_point_location` class implements the incremental randomized algorithm 
introduced by Mulmuley \cgalCite{m-fppa-90} as presented by 
Seidel \cgalCite{s-sfira-91} (see also [\cgalCite{bkos-cgaa-00} Chapter 6). 
It subdivides each arrangement face to pseudo-trapezoidal cells, each 
of constant complexity, and constructs and maintains a linear-size search 
structure on top of these cells, such that each query can be answered 
in \f$ O(\log n)\f$ time, where \f$ n\f$ is the complexity of the arrangement. 

Constructing the search structures takes \f$ O(n \log n)\f$ expected time 
and may require a small number of rebuilds \cgalCite{hkh-iiplgtds-12}. Therefore 
attaching a trapezoidal point-location object to an existing arrangement 
may incur some overhead in running times. In addition, the point-location 
object needs to keep its auxiliary data structures up-to-date as the 
arrangement goes through structural changes. It is therefore recommended 
to use this point-location strategy for static arrangements (or arrangement 
that do not alter frequently), and when the number of issued queries 
is relatively large. 

This strategy supports arbitrary subdivisions, including unbounded ones. 

\cgalModels `ArrangementPointLocation_2`
\cgalModels `ArrangementVerticalRayShoot_2`

\sa `ArrangementPointLocation_2`
\sa `ArrangementVerticalRayShoot_2`
\sa `CGAL::Arr_point_location_result<Arrangement>`
\sa `CGAL_ARR_POINT_LOCATION_VERSION`

*/
template< typename Arrangement >
class Arr_trapezoid_ric_point_location {
public:

/// \name Creation 
/// @{

/*!
If with_guarantees is set to true, the construction performs rebuilds in order to guarantee a resulting structure with linear size and logarithmic query time. Otherwise the structure has expected linear size and expected logarithmic query time. 
*/ 
Arr_trapezoid_ric_point_location (bool with_guarantees = true); 

/*!
Constructs a point location search structure for the given arrangement. If with_guarantees is set to true, the construction performs rebuilds in order to guarantee a resulting structure with linear size and logarithmic query time. Otherwise the structure has expected linear size and expected logarithmic query time. 
*/ 
Arr_trapezoid_ric_point_location (const Arrangement& arr, bool with_guarantees = true); 

/// @} 

/// \name Modifiers 
/// @{

/*!
If with_guarantees is set to true, the structure will guarantee linear size and logarithmic query time, that is, this function may cause a reconstruction of the data structure. 
*/ 
void with_guarantees (bool with_guarantees); 

/// @}

}; /* end Arr_trapezoid_ric_point_location */
} /* end namespace CGAL */
