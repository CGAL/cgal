
/*!
\ingroup PkgArrangementOnSurface2ConceptsDCEL
\cgalConcept

An isolated vertex-information record in a \dcel data structure, which stores
the face that contains the isolated vertex in its interior, along with an
iterator for the isolated vertex in the isolated vertices' container of this
face.

\sa `ArrangementDcel`
\sa `ArrangementDcelFace`

*/

class ArrangementDcelIsolatedVertex {
public:

/// \name Types
/// @{

/*!
the corresponding \dcel face type.
*/
typedef unspecified_type Face;

/*!

*/
typedef Face::Isolated_vertex_iterator Isolated_vertex_iterator;

/// @}

/// \name Creation
/// @{

/*!
default constructor.
*/
Arr_dcel_isolated_vertex();

/// @}

/// \name Access Functions
/// All functions below also have `const` counterparts, returning
/// non-mutable pointers or iterators:
/// @{

/*!
returns the incident face, which contains `iv` in its interior.
*/
Face* face ();

/*!
returns an iterator for the isolated vertex.
*/
Isolated_vertex_iterator iterator();

/// @}

/// \name Modifiers
/// @{

/*!
sets the incident face.
*/
void set_face (Face* f);

/*!
sets the isolated vertex iterator.
*/
void set_iterator (Isolated_vertex_iterator it);

/// @}

}; /* end ArrangementDcelIsolatedVertex */

