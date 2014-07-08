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
}
