// Copyright (c) 1997-2013 INRIA Sophia-Antipolis (France).
// All rights reserved.

namespace CGAL {

/*!
\ingroup PkgPeriodic2Triangulation2VertexFaceClasses

The class `Periodic_2_triangulation_hierarchy_vertex_base_2` is
designed to be used as a vertex base class of a triangulation plugged
into a `Periodic_2_triangulation_hierarchy_2<Tr>`.

It is a model of the concept `TriangulationHierarchyVertexBase_2`
which refines the concept `TriangulationVertexBase_2` and of the
concept `Periodic_2_triangulation_vertex_base_2`.

This class is templated by a parameter `Vb` 
which is to be instantiated by a model of the concept 
`Periodic_2TriangulationVertexBase_2`. 
The class `Periodic_2_triangulation_vertex_base_2<Vb>` inherits 
from the class `Vb`. 
This design allows to use either the default 
vertex base class or a user customized 
vertex base with additional functionalities. 

\cgalModels `TriangulationHierarchyVertexBase_2`
\cgalModels `Periodic_2TriangulationVertexBase_2`

\sa `Periodic_2TriangulationVertexBase_2` 
\sa `TriangulationHierarchyVertexBase_2` 
\sa `CGAL::Periodic_2_triangulation_vertex_base_2<Traits>` 

*/
template< typename Vb >
class Periodic_2_triangulation_hierarchy_vertex_base_2 : public Vb {
public:

/// @}

}; /* end Triangulation_hierarchy_vertex_base_2 */
} /* end namespace CGAL */
