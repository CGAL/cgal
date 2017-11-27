/// Boost Namespace
namespace boost {

/// \ingroup PkgBGLProperties
/// @{

/// The constant `edge_index` is a property tag which identifies the <i>index</i> property of an edge of a \sc{Bgl}
/// <a href="http://www.boost.org/libs/graph/doc/Graph.html"><code>Graph</code></a>.
/// \cgalModels <a href="http://www.boost.org/libs/graph/doc/PropertyTag.html"><code>PropertyTag</code></a>
enum edge_index_t { edge_index};

/// The constant `vertex_index` is a property tag which identifies the <i>index</i> property of a vertex of a \sc{Bgl}
/// <a href="http://www.boost.org/libs/graph/doc/Graph.html"><code>Graph</code></a>.
/// \cgalModels <a href="http://www.boost.org/libs/graph/doc/PropertyTag.html"><code>PropertyTag</code></a>
enum vertex_index_t { vertex_index };

/// The constant `halfedge_index` is a property tag which identifies the <i>index</i> property of a halfedge of a `HalfedgeGraph`.
///
/// This is a property tag introduced by \cgal.
/// \cgalModels <a href="http://www.boost.org/libs/graph/doc/PropertyTag.html"><code>PropertyTag</code></a>
enum halfedge_index_t { halfedge_index };

/// The constant `face_index` is a property tag which identifies the <i>index</i> property of a face of a `FaceGraph`.
///
/// This is a property tag introduced by \cgal.
/// \cgalModels <a href="http://www.boost.org/libs/graph/doc/PropertyTag.html"><code>PropertyTag</code></a>
enum face_index_t { face_index };


/// The constant `vertex_point` is a property tag which refers to the  geometric embedding property of
/// a vertex of a `HalfedgeGraph`.
///
/// This is a property tag introduced by \cgal.
/// \cgalModels <a href="http://www.boost.org/libs/graph/doc/PropertyTag.html"><code>PropertyTag</code></a>
enum vertex_point_t { vertex_point };


/// @}
} // namespace boost

namespace CGAL {


/// @{

  /// \ingroup PkgBGLPropertiesDynamic
/// Dynamic vertex property tag
/// \tparam T the value type of the vertex property
template <typename T>
struct vertex_property_t
{
  /// \param s the name of the property
  vertex_property_t(const std::string s);
};

  /// \ingroup PkgBGLPropertiesDynamic
/// Dynamic halfedge property tag
/// \tparam T the value type of the halfedge property
template <typename T>
struct halfedge_property_t
{
  /// \param s the name of the property
  halfedge_property_t(const std::string s);
};
  
/// \ingroup PkgBGLPropertiesDynamic
/// Helper class to obtain the type of the dynamic property map for a graph type and a dynamic property tag
/// \tparam G must be a model of
/// \tparam D must be a model of DynamicPropertyTag
template <typename G, typename D>
struct dynamic_property_map<G, D>
{
  /// The type of the dynamic property map
  typedef unspecified_type type;
};

/// \ingroup PkgBGLPropertiesDynamic
/// adds a dynamic property map to `g`
template <typename T, typename G>
typename dynamic_property_map<G,vertex_property_t<T> >::const_type
add_property(vertex_property_t<T> prop, const G& g);

  
/// \ingroup PkgBGLPropertiesDynamic
/// removes a dynamic property map assocated to `g`
template<typename G, typename D>
void remove_property(D dpm,
                     const G&);

/// @}
} // namespace CGAL
