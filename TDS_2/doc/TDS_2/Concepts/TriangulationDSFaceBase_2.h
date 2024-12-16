
/*!
\ingroup PkgTDS2Concepts
\cgalConcept

\cgalRefines{TriangulationDataStructure_2::Face}

The concept `TriangulationDSFaceBase_2` describes the requirements for
the face base class of a `CGAL::Triangulation_data_structure_2<Vb,Fb>`.

Note that if the `CGAL::Triangulation_data_structure_2`
is plugged into a triangulation class,
the face base class may have additional geometric
requirements depending on the triangulation class.

At the base level,
(see Sections \ref Section_2D_Triangulations_Software_Design
and \ref TDS_2D_default ),
a face stores handles
on its three vertices and on the three neighboring faces.
The vertices and neighbors are indexed 0,1 and 2.
Neighbor `i` lies opposite to vertex `i`.

Since the `CGAL::Triangulation_data_structure_2` is the class
which defines the handle
types, the face base class has to be somehow
parameterized by the triangulation
data structure. But since the `CGAL::Triangulation_data_structure_2`
itself is parameterized by the face and vertex
base classes, there is a cycle in the definition of these classes.
In order
to break the cycle, the base classes for faces and vertices
which are plugged in to instantiate a
`CGAL::Triangulation_data_structure_2`
use  `void` as triangulation
data structure parameter. Then,
the `CGAL::Triangulation_data_structure_2`
uses a <I>rebind</I> mechanism (similar to the one specified in
`std::allocator`) in order to plug itself
as parameter in the face and vertex base classes.
This mechanism requires that the base class provides
a templated nested class `Rebind_TDS` that
itself provides
the subtype `Rebind_TDS::Other`
which is the <I>rebound</I> version of the base class.
This <I>rebound</I> base class is the class
that the `CGAL::Triangulation_data_structure_2`
actually uses as a base class for the class
`CGAL::Triangulation_data_structure_2::Face`.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Triangulation_ds_face_base_2<TDS>}
\cgalHasModelsEnd

\sa `TriangulationDSVertexBase_2`
\sa `CGAL::Triangulation_data_structure_2<Vb,Fb>`

*/

class TriangulationDSFaceBase_2
{
public:

/// \name Types
/// The concept `TriangulationDSFaceBase_2` has to provide the
/// following types.
/// @{

/*!
This template class has to define a type `Rebind_TDS<TDS_2>::%Other` which is the
<I>rebound</I> face base, where the
`CGAL::Triangulation_data_structure_2` is actually plugged in.
This type `Other` will be the actual base
of the class `CGAL::Triangulation_data_structure_2::Face`.
\note It can be implemented using a nested template class.
\sa Section \ref TDS_2TheRebindMechanism
*/
template <typename TDS2>
using Rebind_TDS = unspecified_type;

/*!

*/
typedef TriangulationDataStructure_2 Triangulation_data_structure;

/*!

*/
typedef TriangulationDataStructure_2::Vertex_handle Vertex_handle;

/*!

*/
typedef TriangulationDataStructure_2::Face_handle Face_handle;

/*!

*/
typedef TriangulationDataStructure_2::Face_data TDS_data;

/// @}

/// \name Creation
/// @{

/*!
default constructor.
*/
TriangulationDSFaceBase_2();

/*!
Initializes the vertices with `v0, v1, v2` and the neighbors
with `Face_handle()`.
*/
TriangulationDSFaceBase_2(Vertex_handle v0,
                          Vertex_handle v1,
                          Vertex_handle v2);

/*!
initializes the vertices with `v0,v1, v2` and the neighbors with
`n0, n1, n2`.
*/
TriangulationDSFaceBase_2(Vertex_handle v0,
                          Vertex_handle v1,
                          Vertex_handle v2,
                          Face_handle n0,
                          Face_handle n1,
                          Face_handle n2);

/// @}

/// \name Access Functions
/// @{

/*!
returns the dimension.
*/
int dimension();

/// @}

/// \name Orientation
/// @{

/*!
Changes the orientation of the face by exchanging `vertex(0)`
with `vertex(1)` and `neighbor(0)` with `neighbor(1)`.
*/
void reorient();

/*!
performs a counterclockwise permutation of the
vertices and neighbors of the face.
*/
void ccw_permute();

/*!
performs a clockwise permutation of the
vertices and neighbors of the face.
*/
void cw_permute();

/// @}

/// \name Checking
/// @{

/*!
performs any required test on a face.

If `verbose` is set to `true`, messages are printed to give
a precise indication of the kind of invalidity encountered.
*/
bool is_valid(bool verbose = false) const;

/// @}

/// \name Various
/// These member functions are required by
/// `CGAL::Triangulation_data_structure_2` because it uses
/// `CGAL::Compact_container` to store its faces. See the documentation of
/// `CGAL::Compact_container` for the exact requirements.
/// @{

/*!

*/
void * for_compact_container() const;

/*!

*/
void for_compact_container(void *p);

/// @}

/// \name Internal
/// \cgalAdvancedBegin
/// These functions are used internally by the triangulation data
/// structure. The user is not encouraged to use them directly as they
/// may change in the future.
/// \cgalAdvancedEnd
/// @{

/*!

*/
TDS_data& tds_data();

/*!

*/
const TDS_data& tds_data() const;

/// @}

}; /* end TriangulationDSFaceBase_2 */

