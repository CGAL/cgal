namespace CGAL {

/*!
\ingroup PkgBGLHelper

The class `Triangulation_face_base_with_id_2` is a model of the
concept `TriangulationFaceBase_2`, the base face of a
2D-triangulation.  It provides an integer field that can be used to
index faces for \sc{Bgl} algorithms.

Note that the user is in charge of setting indices correctly before
running a graph algorithm.

\tparam TriangulationTraits_2 is the geometric traits class
and must be a model of `TriangulationTraits_2`.

\tparam TriangulationFaceBase_2 must be a face base class from which
`Triangulation_face_base_with_id_2` derives. It has the default
value `Triangulation_face_base_2<TriangulationTraits_2>`.

\cgalModels `TriangulationFaceBase_2`

\sa `CGAL::Triangulation_face_base_2`
*/
template< typename TriangulationTraits_2, typename TriangulationFaceBase_2 >
class Triangulation_face_base_with_id_2 : public TriangulationFaceBase_2 {
public:

/// \name Access Functions 
/// @{

/*!
Returns the index. 
*/ 
int id() const; 

/*!
Returns a reference to the index stored in the face.
*/ 
int& id(); 

/// @}

}; /* end Triangulation_face_base_with_id_2 */
} /* end namespace CGAL */
