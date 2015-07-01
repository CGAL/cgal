
namespace CGAL {

/*!
\ingroup PkgMesh2



The class `Delaunay_mesh_face_base_2` is a model for the concept 
`DelaunayMeshFaceBase_2`. 

This class can be used directly or it can serve as a base to derive other 
classes with some additional attributes (a color for example) tuned to a 
specific application. 



\tparam Traits is the geometric traits class. It must 
be the same as the one used for the Delaunay mesh. 

\tparam Fb is the base class from which `Delaunay_mesh_face_base_2` 
derives. It must be a model of the `ConstrainedTriangulationFaceBase_2` concept. 


\cgalModels `DelaunayMeshFaceBase_2`

*/
template< typename Traits, typename Fb >
class Delaunay_mesh_face_base_2 : Fb {
public:

/// @}

}; /* end Delaunay_mesh_face_base_2 */
} /* end namespace CGAL */
