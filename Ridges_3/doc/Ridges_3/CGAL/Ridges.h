namespace CGAL {

/*!
\ingroup PkgRidges3Enums
The enum `Ridge_order` describes the order of differential quantities
used to distinguish elliptic and hyperbolic ridges. Third or fourth
order quantities may be used as explained in Section  \ref Ridges_3Approximating
of the User Manual.

\sa `Ridge_approximation`
*/
enum Ridge_order { Ridge_order_3 = 3, Ridge_order_4};

/*!
\ingroup PkgRidges3Enums
The enum `Ridge_type` describes the types for the class `Ridge_line`.
\sa `Ridge_line`
*/
enum Ridge_type {
  MAX_ELLIPTIC_RIDGE = 1, MAX_HYPERBOLIC_RIDGE,
  MAX_CREST_RIDGE, MIN_ELLIPTIC_RIDGE,
  MIN_HYPERBOLIC_RIDGE, MIN_CREST_RIDGE
};


/*!
\ingroup PkgRidges3Ref

The function `compute_crest_ridges()` is a shortcut to the
method `Ridge_approximation::compute_crest_ridges()`.
See `Ridge_approximation::Ridge_approximation()` for an explanation of the parameters.


*/
template < class TriangleMesh,
           class VertexFTMap,
           class VertexVectorMap,
           class OutputIterator>
OutputIterator compute_crest_ridges(const TriangleMesh &tm,
                                    VertexFTMap vertex_k1_pm,
                                    VertexFTMap vertex_k2_pm,
                                    VertexFTMap vertex_b0_pm,
                                    VertexFTMap vertex_b3_pm,
                                    VertexVectorMap vertex_d1_pm,
                                    VertexVectorMap vertex_d2_pm,
                                    VertexFTMap vertex_P1_pm,
                                    VertexFTMap vertex_P2_pm,
                                    OutputIterator it,
                                    CGAL::Ridge_order order = CGAL::Ridge_order_3);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgRidges3Ref

The function `compute_max_ridges()` is a shortcut to the
method `Ridge_approximation::compute_max_ridges()`.
See `Ridge_approximation::Ridge_approximation()` for an explanation of the parameters.
*/
template < class TriangleMesh,
           class VertexFTMap,
           class VertexVectorMap,
           class OutputIterator>
OutputIterator compute_max_ridges(const TriangleMesh &tm,
                                  VertexFTMap vertex_k1_pm,
                                  VertexFTMap vertex_k2_pm,
                                  VertexFTMap vertex_b0_pm,
                                  VertexFTMap vertex_b3_pm,
                                  VertexVectorMap vertex_d1_pm,
                                  VertexVectorMap vertex_d2_pm,
                                  VertexFTMap vertex_P1_pm, const
                                  VertexFTMap vertex_P2_pm,
                                  OutputIterator it,
                                  CGAL::Ridge_order order = CGAL::Ridge_order_3);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgRidges3Ref

The function `compute_min_ridges()` is a shortcut to
the method `Ridge_approximation::compute_min_ridges()`.
See `Ridge_approximation::Ridge_approximation()` for an explanation of the parameters.


*/
template < class TriangleMesh,
class VertexFTMap,
class VertexVectorMap,
class OutputIterator>
OutputIterator compute_min_ridges(const TriangleMesh &tm,
                                  VertexFTMap vertex_k1_pm,
                                  VertexFTMap vertex_k2_pm,
                                  VertexFTMap vertex_b0_pm,
                                  VertexFTMap vertex_b3_pm,
                                  VertexVectorMap vertex_d1_pm,
                                  VertexVectorMap vertex_d2_pm,
                                  VertexFTMap vertex_P1_pm,
                                  VertexFTMap vertex_P2_pm,
                                  OutputIterator it,
                                  CGAL::Ridge_order order = CGAL::Ridge_order_3);

} /* namespace CGAL */


namespace CGAL {

/*!
\ingroup PkgRidges3Ref

The class `Ridge_approximation` computes the approximation of
ridges of a triangular polyhedral surface.

\tparam TriangleMesh is the surface type. In the following let `K` be `Kernel_traits<boost::property_traits<TriangleMesh,CGAL::vertex_point_t>::%value_type>::%Kernel`
\tparam VertexFTMap A property map with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and `K::FT` as value type.
\tparam VertexVectorMap A property map with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and `K::Vector_3` as value type.

\pre (checked at compile time)
\pre The types `K::FT` and
 `boost::property_traits<VertexFTMap>::%value_type` must coincide.
\pre The types `K::Vector_3` and
  `boost::property_traits<VertexVectorMap>::%value_type` must coincide.
\pre The types `boost::graph_traits<TriangleMesh>::%vertex_descriptor`,
  and`boost::property_traits<VertexFTMap>::%key_type`, and
  `boost::property_traits<VertexVectorMap>::%key_type` must coincide.

\sa `Ridge_line`

*/
template< typename TriangleMesh, typename VertexFTMap, typename VertexVectorMap >
class Ridge_approximation {
public:


/// \name Creation
/// @{

/*!
The two last property maps may
not be used if computations are performed with
the parameter `Ridges_order_3`, in which case these
property maps shall be initialized with their
default constructors.

\param tm the triangle mesh
\param vertex_k1_pm maximal principal curvatures
\param vertex_k2_pm minimal principal curvatures
\param vertex_b0_pm third order extremalities
\param vertex_b3_pm third order extremalities
\param vertex_d1_pm maximal principal directions of curvature
\param vertex_d2_pm minimal principal directions of curvature
\param vertex_P1_pm fourth order quantities
\param vertex_P2_pm fourth order quantities
*/
Ridge_approximation(const TriangleMesh &tm,
                    VertexFTMap vertex_k1_pm,
                    VertexFTMap vertex_k2_pm,
                    VertexFTMap vertex_b0_pm,
                    VertexFTMap vertex_b3_pm,
                    VertexVectorMap vertex_d1_pm,
                    VertexVectorMap vertex_d2_pm,
                    VertexFTMap vertex_P1_pm,
                    VertexFTMap vertex_P2_pm);

/*!
Outputs ridges of types `MAX_ELLIPTIC_RIDGE` and `MAX_HYPERBOLIC_RIDGE`.
\tparam OutputIterator an output iterator wìth value type `Ridge_line*`.
*/
template <class OutputIterator> OutputIterator compute_max_ridges(OutputIterator it, Ridge_order ord = Ridge_order_3);

/*!
Outputs ridges of types `MIN_ELLIPTIC_RIDGE` and `MIN_HYPERBOLIC_RIDGE`.
\tparam OutputIterator an output iterator with
value type `Ridge_line*`.
*/
template <class OutputIterator> OutputIterator compute_min_ridges(OutputIterator it, Ridge_order ord = Ridge_order_3);

/*!
Outputs ridges of types `MAX_CREST_RIDGE` and `MIN_CREST_RIDGE`.
\tparam OutputIterator is an output iterator with
value type `Ridge_line*`.
*/
template <class OutputIterator> OutputIterator compute_crest_ridges(OutputIterator it, Ridge_order ord = Ridge_order_3);

/// @}

}; /* end Ridge_approximation */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgRidges3Ref

The class `Ridge_line` stores the description of a ridge line. The list
of halfedges defines a connected sequence of edges (but not as
oriented halfedges). The scalar \f$ b\f$ paired with a halfedge \f$ pq\f$ is the
barycentric coordinate of the crossing point \f$ r\f$ with the ridge:
\f$ r = b\times p + (1-b)\times q\f$.

\sa `Ridge_approximation`

*/
template< typename TriangleMesh >
class Ridge_line {
public:

/// \name Types
/// @{

/*!

*/
typedef typename TriangleMesh::Traits::FT FT;

/*!

*/
  typedef typename boot::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

/*!
A halfedge crossed by a ridge is paired with the barycentric
coordinate of the crossing point.
*/
typedef std::pair< halfedge_descriptor, FT> Ridge_halfedge;

/// @}

/// \name Creation
/// @{

/*!
default constructor.
*/
Ridge_line();

/// @}

/// \name Access Functions
/// @{

/*!

*/
Ridge_type line_type() const;

/*!

*/
FT strength() const;

/*!

*/
FT sharpness() const;

/*!

*/
const std::list<Ridge_halfedge>* line() const;


/// @}

}; /* end Ridge_line */

/*!
\relates Ridge_line
Writes the line type, strength, sharpness and coordinates of the
points of the polyline to `os`.
*/
template< typename TriangleMesh >
std::ostream& operator<<(std::ostream& os, const Ridge_line<TriangleMesh>& r);

} /* end namespace CGAL */

