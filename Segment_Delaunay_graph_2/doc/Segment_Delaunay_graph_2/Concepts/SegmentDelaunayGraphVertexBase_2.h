
/*!
\ingroup PkgSegmentDelaunayGraph2Concepts
\cgalConcept

\cgalRefines{TriangulationDSVertexBase_2}

The concept `SegmentDelaunayGraphVertexBase_2` describes the
requirements for the vertex base class of the
`SegmentDelaunayGraphDataStructure_2` concept. A vertex stores a
site of the segment Delaunay graph and provides access to one of its
incident faces through a `Face_handle`.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Segment_Delaunay_graph_vertex_base_2<St>}
\cgalHasModelsEnd

\sa `SegmentDelaunayGraphTraits_2`
\sa `SegmentDelaunayGraphSite_2`
\sa `SegmentDelaunayGraphStorageSite_2`
\sa `CGAL::Segment_Delaunay_graph_site_2<K>`
\sa `CGAL::Segment_Delaunay_graph_storage_site_2<Gt,SSTag>`
\sa `CGAL::Triangulation_data_structure_2<Vb,Fb>`
*/
class SegmentDelaunayGraphVertexBase_2 {
public:

/// \name Types
/// @{

/*!
A type for the geometric traits that defines the site.
\pre The type `Geom_traits` must define the type `Site_2`.
*/
typedef unspecified_type Geom_traits;

/*!
A type for the site. This type must coincide with the type `Geom_traits::Site_2`.
*/
typedef unspecified_type Site_2;

/*!
A type that indicates what kind of storage type to use.
`Storage_site_tag` must either be `CGAL::Tag_true` or `CGAL::Tag_false`.
*/
typedef unspecified_type Storage_site_tag;

/*!
A type for the internal representation of sites.
This type must satisfy the requirements of the concept `SegmentDelaunayGraphStorageSite_2`.
*/
typedef unspecified_type Storage_site_2;

/*!
A type for the underlying data structure, to which the vertex belongs to.
*/
typedef unspecified_type Data_structure;

/*!
A type for the vertex handle of the segment Delaunay graph data structure.
*/
typedef unspecified_type Vertex_handle;

/*!
A type for the face handle of the segment Delaunay graph data structure.
*/
typedef unspecified_type Face_handle;

/// @}

/// \name Creation
/// In addition to the default and copy constructors and following
/// constructors are required:
/// @{

/*!
Constructs a vertex associated with the site represented by the storage site `ss`.
*/
SegmentDelaunayGraphVertexBase_2(const Storage_site_2& ss);

/*!
Constructs a vertex associated with the site represented by the storage site `ss`,
and pointing to the face associated with the face handle `f`.
*/
SegmentDelaunayGraphVertexBase_2(const Storage_site_2& ss, Face_handle f);

/// @}

/// \name Access Functions
/// @{

/*!
Returns the storage site representing the site.
*/
const Storage_site_2& storage_site();

/*!
Returns the site.
*/
Site_2 site();

/*!
Returns a handle to an incident face.
*/
Face_handle face();

/// @}

/// \name Setting
/// @{

/*!
Sets the storage site.
*/
void set_site(const Storage_site_2& ss);

/*!
Sets the incident face.
*/
void set_face(Face_handle f);

/// @}

/// \name Checking
/// @{

/*!
Performs any required tests on a vertex.
*/
bool is_valid(bool verbose, int level) const;

/// @}

}; /* end SegmentDelaunayGraphVertexBase_2 */

