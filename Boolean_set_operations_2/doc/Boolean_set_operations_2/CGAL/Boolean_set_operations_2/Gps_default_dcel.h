
namespace CGAL {

/*!
\ingroup PkgBooleanSetOperations2

The class `Gps_default_dcel` is the default <span class="textsc">Dcel</span> class used by the 
`General_polygon_set_2` and `Polygon_set_2` class-templates 
to represent the underlying internal `Arrangement_2` data structure. 

\cgalModels ::GeneralPolygonSetDcel 

*/
template< typename Traits >
class Gps_default_dcel 
  : CGAL::Arr_dcel_base< Arr_vertex_base<typename Traits_::Point_2>,
                         Arr_halfedge_base<typename Traits_::X_monotone_curve_2>,
                         Gps_face_base> {
public:

/// \name Types 
/// @{

/*! 
allows the rebinding of the <span class="textsc">Dcel</span> with a different traits class `T`. 
*/ 
typedef Hidden_type template <class T> rebind; 

/// @}

}; /* end Gps_default_dcel */
} /* end namespace CGAL */
