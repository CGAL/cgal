/// CGAL Namespace
namespace CGAL {

/// \ingroup PkgBGLProperties
/// @{

/// The constant `vertex_index` is a property tag which identifies the <i>index</i> property of a vertex of a \bgl
/// <a href="https://www.boost.org/libs/graph/doc/Graph.html"><code>Graph</code></a>.
/// \cgalModels <a href="https://www.boost.org/libs/graph/doc/PropertyTag.html"><code>PropertyTag</code></a>
enum vertex_index_t { vertex_index };

/// The constant `halfedge_index` is a property tag which identifies the <i>index</i> property of a halfedge of a `HalfedgeGraph`.
///
/// This is a property tag introduced by \cgal.
/// \cgalModels <a href="https://www.boost.org/libs/graph/doc/PropertyTag.html"><code>PropertyTag</code></a>
enum halfedge_index_t { halfedge_index };

/// The constant `edge_index` is a property tag which identifies the <i>index</i> property of an edge of a \bgl
/// <a href="https://www.boost.org/libs/graph/doc/Graph.html"><code>Graph</code></a>.
/// \cgalModels <a href="https://www.boost.org/libs/graph/doc/PropertyTag.html"><code>PropertyTag</code></a>
enum edge_index_t { edge_index };

/// The constant `face_index` is a property tag which identifies the <i>index</i> property of a face of a `FaceGraph`.
///
/// This is a property tag introduced by \cgal.
/// \cgalModels <a href="https://www.boost.org/libs/graph/doc/PropertyTag.html"><code>PropertyTag</code></a>
enum face_index_t { face_index };

/// The constant `vertex_point` is a property tag which refers to the  geometric embedding property of
/// a vertex of a `HalfedgeGraph`.
///
/// This is a property tag introduced by \cgal.
/// \cgalModels <a href="https://www.boost.org/libs/graph/doc/PropertyTag.html"><code>PropertyTag</code></a>
enum vertex_point_t { vertex_point };

/// @}

/// \ingroup PkgBGLProperties
///
/// \brief graph_has_property is used to indicate if a model of `HalfedgeGraph` or `FaceGraph`
/// has an internal property associated with the given `PropertyTag`.
///
/// It inherits from \link Tag_true `CGAL::Tag_true` \endlink if there is a
/// default internal property map for the corresponding property tag and from
/// \link Tag_false `CGAL::Tag_false` \endlink otherwise.
///
/// \tparam Graph a model of `HalfedgeGraph` or `FaceGraph`
/// \tparam PropertyTag the type of a property tag referring to the property of interest.
///
template<typename Graph, typename PropertyTag>
struct graph_has_property;

/// @{

/// \ingroup PkgBGLPropertiesDynamic
/// Dynamic vertex property tag
/// \tparam T the value type of the vertex property
template <typename T>
struct dynamic_vertex_property_t
{
  dynamic_vertex_property_t();
};

/// \ingroup PkgBGLPropertiesDynamic
/// Dynamic halfedge property tag
/// \tparam T the value type of the halfedge property
template <typename T>
struct dynamic_halfedge_property_t
{
  dynamic_halfedge_property_t();
};

/// \ingroup PkgBGLPropertiesDynamic
/// Dynamic edge property tag
/// \tparam T the value type of the edge property
template <typename T>
struct dynamic_edge_property_t
{
  dynamic_edge_property_t();
};

/// \ingroup PkgBGLPropertiesDynamic
/// Dynamic face property tag
/// \tparam T the value type of the face property
template <typename T>
struct dynamic_face_property_t
{
  dynamic_face_property_t();
};

/// @}
} // namespace CGAL
