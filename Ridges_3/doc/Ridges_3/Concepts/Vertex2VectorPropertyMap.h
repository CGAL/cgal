
/*!
\ingroup PkgRidges_3Concepts
\cgalConcept

The concept `Vertex2VectorPropertyMap` specializes the concept of `LvaluePropertyMap`.
It is intended to be used in combination with 
the concept `TriangulatedSurfaceMesh` in the class 
`CGAL::Ridge_approximation`. It associates a three dimensional vector 
`TriangulatedSurfaceMesh::Traits::Vector_3` to keys which are 
`TriangulatedSurfaceMesh::Vertex_handle`. 

\cgalHasModel `CGAL::Vertex2Data_Property_Map_with_std_map::Vertex2Vector_property_map`

\sa `TriangulatedSurfaceMesh`
*/
class Vertex2VectorPropertyMap {
public:

/// @}

}; /* end Vertex2VectorPropertyMap */

