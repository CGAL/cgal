
namespace CGAL {

/*!
\ingroup PkgArrangement2PointLocation

\anchor arr_reflm_pl 

The `Arr_landmarks_point_location` class implements a Jump & Walk algorithm, where special 
points, referred to as "landmarks", are chosen in a preprocessing stage, 
their place in the arrangement is found, and they are inserted into a 
data-structure that enables efficient nearest-neighbor search (a 
<span class="textsc">Kd</span>-tree). Given a query point, the nearest landmark is located and a 
"walk" strategy is applied from the landmark to the query point. 

There are various strategies to select the landmark set in the 
arrangement, where the strategy is determined by the 
`Generator` template parameter. The following landmark-generator 
classes are available: 
<DL> 
<DT><B>`Arr_landmarks_vertices_generator` - </B><DD> 
The arrangement vertices are used as the landmarks set. 

<DT><B>`Arr_random_landmarks_generator` - </B><DD> 
\f$ n\f$ random points in the bounding box of the arrangement are chosen 
as the landmarks set. 

<DT><B>`Arr_halton_landmarks_generator` - </B><DD> 
\f$ n\f$ Halton points in the bounding box of the arrangement are chosen 
as the landmarks set. 

<DT><B>`Arr_middle_edges_landmarks_generator` - </B><DD> 
The midpoint of each arrangement edge is computed, and the resulting 
set of points is used as the landmarks set. This generator can be applied 
only for arrangements of line segments. 

<DT><B>`Arr_grid_landmarks_generator` - </B><DD> 
A set of \f$ n\f$ landmarks are chosen on a 
\f$ \lceil \sqrt{n} \rceil \times \lceil \sqrt{n} \rceil\f$ 
grid that covers the bounding box of the arrangement. 
</DL> 
The `Arr_landmarks_vertices_generator` class is the default generator 
and used if no `Generator` parameter is specified. 

It is recommended to use the landmarks point-location strategy 
when the application frequently issues point-location queries on a 
rather static arrangement that the changes applied to it are mainly 
insertions of curves and not deletions of them. 

\cgalModels `ArrangementPointLocation_2`
\cgalModels `ArrangementVerticalRayShoot_2`

\sa `ArrangementPointLocation_2`
\sa `ArrangementVerticalRayShoot_2`
\sa `CGAL::Arr_point_location_result<Arrangement>`
\sa `CGAL_ARR_POINT_LOCATION_VERSION`

*/
template< typename Arrangement, typename Generator >
class Arr_landmarks_point_location {
public:

/// @}

}; /* end Arr_landmarks_point_location */
} /* end namespace CGAL */
