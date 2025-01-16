
/*!
\ingroup PkgApolloniusGraph2Concepts
\cgalConcept

The concept `ApolloniusGraphVertexBase_2` describes the
requirements for the vertex base class of the
`ApolloniusGraphDataStructure_2` concept. A vertex stores an
Apollonius site and provides access to one of its incident faces
through a `Face_handle`. In addition, it maintains a container of
sites. The container stores the hidden sites related to the vertex.

\cgalRefines{TriangulationVertexBase_2}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Apollonius_graph_vertex_base_2<Gt,StoreHidden>}
\cgalHasModelsEnd

\sa `ApolloniusGraphDataStructure_2`
\sa `CGAL::Apollonius_graph_2<Gt,Agds>`
\sa `CGAL::Triangulation_data_structure_2<Vb,Fb>`
*/
class ApolloniusGraphVertexBase_2 {
public:

/// \name Types
/// @{

/*!
A type for the geometric traits that defines
the site stored. \pre The type `Geom_traits` must define the type `Site_2`.
*/
typedef unspecified_type Geom_traits;

/*!
A Boolean that indicates if
hidden sites are actually stored or not. Its value is `true` if
hidden sites are stored, `false` otherwise.
*/
typedef unspecified_type Store_hidden;

/*!
A type for the site stored.
\pre This type must coincide with the type `Geom_traits::Site_2`.
*/
typedef unspecified_type Site_2;

/*!
A type for the
Apollonius graph data structure, to which the vertex belongs to.
*/
typedef unspecified_type Apollonius_graph_data_structure_2;

/*!
A type for the vertex handle of the
Apollonius graph data structure.
*/
typedef unspecified_type Vertex_handle;

/*!
A type for the face handle of the
Apollonius graph data structure.
*/
typedef unspecified_type Face_handle;

/*!
An iterator that
iterates over the hidden sites in the hidden sites
container of the vertex.
\pre Must be a model of `Iterator`.
*/
typedef unspecified_type Hidden_sites_iterator;

/// @}

/// \name Creation
/// @{

/*!
%Default constructor.
*/
ApolloniusGraphVertexBase_2();

/*!
Constructs a vertex associated with the Apollonius site `s` and
embedded at the center of `s`.
*/
ApolloniusGraphVertexBase_2(Site_2 s);

/*!
Constructs a vertex associated with
the site `s`, embedded at the center of `s`,
and pointing to the face associated with the face handle `f`.
*/
ApolloniusGraphVertexBase_2(Site_2 s,
Face_handle f);

/// @}

/// \name Access Functions
/// @{

/*!
Returns the Apollonius site.
*/
Site_2 site();

/*!
Returns a handle to an incident face.
*/
Face_handle face();

/*!
Returns the number of hidden sites in the hidden
sites container.
*/
unsigned int number_of_hidden_sites();

/*!
Starts at an arbitrary hidden site.
*/
Hidden_sites_iterator
hidden_sites_begin();

/*!
Past-the-end iterator.
*/
Hidden_sites_iterator hidden_sites_end();

/// @}

/// \name Setting and unsetting
/// @{

/*!
Sets the Apollonius site.
*/
void set_site(Site_2 s);

/*!
Sets the incident face.
*/
void set_face(Face_handle f);

/*!
Adds a hidden site to the container of hidden sites.
*/
void add_hidden_site(Site_2 s);

/*!
Clears the container of hidden sites.
*/
void clear_hidden_sites_container();

/// @}

/// \name Checking
/// @{

/*!
Performs any required tests on a vertex.
*/
bool is_valid(bool verbose, int level) const;

/// @}

}; /* end ApolloniusGraphVertexBase_2 */

