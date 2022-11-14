
namespace CGAL {

/*!
\ingroup PkgApolloniusGraph2Ref

The class `Apollonius_graph_vertex_base_2` provides a model for the
`ApolloniusGraphVertexBase_2` concept which is the vertex base
required by the `ApolloniusGraphDataStructure_2` concept. The
class `Apollonius_graph_vertex_base_2` has two template arguments, the first being the
geometric traits of the Apollonius graph and should be a model of the
concept `ApolloniusGraphTraits_2`. The second is a Boolean which
controls whether hidden sites are actually stored. Such a
control is important if the user is not interested in hidden sites
and/or if only insertions are made, in which case no hidden
site can become visible. If `StoreHidden` is set to
`true`, hidden sites are stored, otherwise they are
discarded. By default `StoreHidden` is set to `true`.

\cgalModels `ApolloniusGraphVertexBase_2`

\sa `CGAL::Triangulation_data_structure_2<Vb,Fb>`
\sa `CGAL::Apollonius_graph_hierarchy_vertex_base_2<Gt>`
*/
template< typename Gt, typename StoreHidden >
class Apollonius_graph_vertex_base_2 {
public:

/// \name Creation
/// @{

/*!
%Default constructor.
*/
Apollonius_graph_vertex_base_2();

/*!
Constructs a vertex associated with the site `s` and
embedded at the center of `s`.
*/
Apollonius_graph_vertex_base_2(const Site_2& s);

/*!
Constructs a vertex associated with
the site `s`, embedded at the center of `s`,
and pointing to the face associated with the face
handle `f`.
*/
Apollonius_graph_vertex_base_2(const Site_2& s, Face_handle f);

/// @}

}; /* end Apollonius_graph_vertex_base_2 */
} /* end namespace CGAL */
