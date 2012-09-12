
namespace CGAL {

/*!
\ingroup PkgArrangement2

\anchor arr_refnaive_pl 

The `Arr_naive_point_location` class implements a naive algorithm that traverses 
all the vertices and halfedges in the arrangement in search for an 
answer to a point-location query. 
The query time is therefore linear in the complexity of the arrangement. 
Naturally, this point-location strategy could turn into a heavy 
time-consuming process when applied to dense arrangements. 

\models ::ArrangementPointLocation_2 
\models ::ArrangementVerticalRayShoot_2 

*/
template< typename Arrangement >
class Arr_naive_point_location {
public:

/// @}

}; /* end Arr_naive_point_location */
} /* end namespace CGAL */
