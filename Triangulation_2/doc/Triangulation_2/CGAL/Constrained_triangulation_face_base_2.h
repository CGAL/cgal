
namespace CGAL {

/*!
\ingroup PkgTriangulation2VertexFaceClasses

The class `Constrained_triangulation_face_base_2` is the default model for the concept 
`ConstrainedTriangulationFaceBase_2` to be used as base face class 
of constrained triangulations. 

\models ::ConstrainedTriangulationFaceBase_2 

Parameters 
-------------- 

The first template parameter is a geometric traits. 

The second template parameter has to be a model 
of the concept `TriangulationFaceBase_2`. 
Its default is `CGAL::Triangulation_face_base_2<Traits>` 

Inherits From 
-------------- 

The class `Constrained_triangulation_face_base_2` derives from its
parameter `Fb`.  and add three Boolean to deal with information about
constrained edges.

The member functions `cw(int i)`, `ccw(int i)` 
and `reorient` are overloaded to update 
information about constrained edges. 

\sa `TriangulationFaceBase_2` 
\sa `ConstrainedTriangulationFaceBase_2` 
\sa `CGAL::Constrained_triangulation_2<Traits,Tds>` 
\sa `CGAL::Triangulation_face_base_2<Traits>` 

*/
template< typename Traits, typename Fb >
class Constrained_triangulation_face_base_2 : public Fb {
public:

/// @}

}; /* end Constrained_triangulation_face_base_2 */
} /* end namespace CGAL */
