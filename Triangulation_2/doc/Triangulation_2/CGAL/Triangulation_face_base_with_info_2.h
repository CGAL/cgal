
namespace CGAL {

/*!
\ingroup PkgTriangulation2VertexFaceClasses

The class `Triangulation_face_base_with_info_2` is a model of the concept 
`TriangulationFaceBase_2` to be plugged into the 
triangulation data structure of a triangulation class. 
It provides an easy way to add some user defined information 
in the faces of a triangulation. 

### Parameters ###

The first template argument is the information the user would like to add 
to a face. It has to be `DefaultConstructible` and `Assignable`. 

The second template argument is a geometric traits class 
and is actually not used in `Triangulation_face_base_with_info_2` . 

The third parameter is a face base class from which 
`Triangulation_face_base_with_info_2` derives. 

\models Because `Triangulation_face_base_with_info_2` derives from the class instantiating its third 
parameter, it will be a model of the same face base concept as its parameter: 
`TriangulationFaceBase_2`,  `ConstrainedTriangulationFaceBase_2`, or ::RegularTriangulationFaceBase_2 

\sa `CGAL::Triangulation_face_base_2<Traits,Fb>` 
\sa `CGAL::Constrained_triangulation_face_base_2<Traits,Fb>` 
\sa `CGAL::Regular_triangulation_face_base_2<Traits,Fb>` 

*/
template< typename Info, typename Traits, typename Fb >
class Triangulation_face_base_with_info_2 : public Fb {
public:

/// \name Types 
/// @{

/*! 

*/ 
typedef Info Info; 

/// @} 

/// \name Access Functions 
/// @{

/*! 
Returns a const reference to the object of type `Info` stored in the face. 
*/ 
const Info& info() const; 

/*! 
Returns a reference to the object of type `Info` stored in the face. 
*/ 
Info & info(); 

/// @}

}; /* end Triangulation_face_base_with_info_2 */
} /* end namespace CGAL */
