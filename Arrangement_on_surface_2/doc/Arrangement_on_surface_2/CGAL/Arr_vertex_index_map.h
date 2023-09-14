namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2Ref

`Arr_vertex_index_map` maintains a mapping of vertex handles of an
attached arrangement object to indices (of type `unsigned int`).
This class template is a model of the concept
`ReadablePropertyMap`. A mapping between vertex handles and indices
enables convenient usage of property-map classes supplied by `boost`.
For example, the property-map class templates
`boost::vector_property_map`, which is based on `std::vector`,
and `boost::iterator_property_map`, which can be used to implement
a property map based on a native \CC array, require the
user to supply a mapping such as `Arr_vertex_index_map`.

As new vertices might be inserted into the attached arrangement, and
existing vertices might be removed, the notification mechanism is used
to dynamically maintain the mapping of vertex handles to indices.


\cgalModels{DefaultConstructible,CopyConstructible,Assignable,ReadablePropertyMap}

\sa `Arr_observer<Arrangement>`
\sa `Arr_face_index_map<Arrangement>`
*/

template< typename Arrangement >
class Arr_vertex_index_map: public Arr_observer<Arrangement> {
public:

/// \name Types
/// @{

/*!
the type of the attached arrangement.
*/
typedef Arrangement Arrangement_2;

typedef boost::readable_property_map_tag category;

typedef unsigned int value_type;

typedef unsigned int reference;

typedef Vertex_handle key_type;

/*!
The vertex handle type.
*/
typedef typename Arrangement_2::Vertex_handle Vertex_handle;

/*!
The type of mapping of vertices to indices.
*/
typedef Unique_hash_map<Vertex_handle, value_type> Index_map;

/// @}

/// \name Creation
/// @{

/*!
constructs a map that is unattached to any arrangement instance.
*/
Arr_vertex_index_map();

/*!
constructs a map and attaches it to the given arrangement `arr`.
*/
Arr_vertex_index_map(Arrangement_2& arr);

/// @}

}; /* end Arr_accessor */
} /* end namespace CGAL */






