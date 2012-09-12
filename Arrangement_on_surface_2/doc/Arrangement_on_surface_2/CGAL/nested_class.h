
/*!


/*!


The basic <span class="textsc">Dcel</span> face type. Serves as a basis class for an extended 
face record with auxiliary data fields. 

\models ::ArrangementDcelFace 

*/
class Arr_dcel_base::Arr_face_base {

}; /* end Arr_dcel_base::Arr_face_base */

/*!


The basic <span class="textsc">Dcel</span> halfedge type. Serves as a basis class for an 
extended halfedge record with auxiliary data fields. The `Curve` 
parameter is the type of \f$ x\f$-monotone curves associated with the vertices. 

\models ::ArrangementDcelHalfedge 

*/
template< typename Curve >
class Arr_dcel_base::Arr_halfedge_base {

}; /* end Arr_dcel_base::Arr_halfedge_base */

/*!


The basic <span class="textsc">Dcel</span> vertex type. Serves as a basis class for an extended 
vertex record with auxiliary data fields. The `Point` parameter is 
the type of points associated with the vertices. 

\models ::ArrangementDcelVertex 

*/
template< typename Point >
class Arr_dcel_base::Arr_vertex_base {

}; /* end Arr_dcel_base::Arr_vertex_base */

