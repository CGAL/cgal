
/*!
\ingroup PkgRidges_3Concepts
\cgalConcept

The concept `Vertex2FTPropertyMap` specializes the concept of `LvaluePropertyMap`.
It is intended to be used in combination with 
the concept `TriangulatedSurfaceMesh` in the class 
`CGAL::Ridge_approximation`. It associates a field type value 
`TriangulatedSurfaceMesh::Traits::FT` to keys which are 
`TriangulatedSurfaceMesh::Vertex_handle`. 

\cgalHasModel `CGAL::Vertex2Data_Property_Map_with_std_map::Vertex2FT_property_map`

\sa `TriangulatedSurfaceMesh`

*/

class Vertex2FTPropertyMap {
public:

/// @}

}; /* end Vertex2FTPropertyMap */

