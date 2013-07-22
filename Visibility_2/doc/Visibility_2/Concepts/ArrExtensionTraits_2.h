/*!
\ingroup PkgVisibility_2Concepts
\cgalConcept

All visibility polgyon algorithms provided in \cgal are parameterized with a traits class 'Traits', which defines the extension of Arrangement_2 the output will have.

\cgalHasModel `CGAL::Arr_extension_default_traits_2<K>`

\sa `Visibility_2`

*/
class ArrExtensionTraits_2 {
public :



/// \name Types
/// @{

/*!
 * The vertex handle type on which the functors will operate.
 */
typedef Hidden_type Vertex_handle;

/*!
 * The halfedge handle type on which the functors will operate.
 */
typedef Hidden_type Halfedge_handle;

/*!
 * The face handle type on which the functors will operate.
 */
typedef Hidden_type Face_handle;

/*!
 * Add auxiliary information to vertex.
 */
typedef Hidden_type Extend_vertex;

/*!
 * Add auxiliary information to halfedge.
 */
typedef Hidden_type Extend_halfedge;

/*!
 * Add auxiliary information to face.
 */
typedef Hidden_type Extend_face;

/// \name Creation
/// @{
/*!
default creator
*/
ArrExtensionTraits_2 ();

/*!
copy creator
*/
ArrExtensionTraits_2 (const ArrExtensionTraits_2 & Traits);

/// @}

/// \name Operations
/// The following member functions to create instances of the above predicate oject types.
/// @{

/*!
 *
 */
Extend_vertex extend_vertex_object();

/*!
 *
 */
Extend_halfedge extend_halfedge_object();

/*!
 *
 */
Extend_face extend_face_object();

/// @}


}

