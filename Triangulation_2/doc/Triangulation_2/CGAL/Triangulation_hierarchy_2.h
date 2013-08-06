
namespace CGAL {

/*!
\ingroup PkgTriangulation2TriangulationClasses

The class `Triangulation_hierarchy_2` implements a triangulation augmented with 
a data structure which allows fast point location queries. 

The data structure is a hierarchy 
of triangulations. The triangulation at the lowest level is 
the original triangulation where operations and point location are to 
be performed. 
Then at each succeeding level, the data structure 
stores a triangulation of a small random sample of the vertices 
of the triangulation at the preceding level. 

Point location 
is done through a top-down nearest neighbor query. 
The nearest neighbor query is first 
performed naively in the top level triangulation. 
Then, at each following level, the nearest neighbor at that level 
is found through a linear walk performed from 
the nearest neighbor found at the preceding level. 

Because the number of vertices in each triangulation is only a small 
fraction of the number of vertices of the preceding triangulation 
the data structure remains small and achieves fast point location 
queries on real 
data. As proved in \cgalCite{d-iirdt-98}, this structure has an optimal behavior 
when it is built for Delaunay triangulations. 
However it can be used as well for other triangulations. 


\tparam Tr may be any of the \cgal triangulation classes. 

\cgalHeading{Types}

The class `Triangulation_hierarchy_2` inherits the types from its base triangulation 
class `Tr`. 

The class `Triangulation_hierarchy_2` offers exactly the same functionalities 
as the triangulation Tr does. 
Location queries are overloaded to benefit from the 
data structure. Modifiers (insertion, removal, and displacement) are overloaded 
to take care of updating the data structure. 

Be careful that I/O operations are not overloaded. 
Writing a `Triangulation_hierarchy_2` into a file 
writes only the lowest level triangulation and drops the hierarchy 
and reading it from a file results in a triangulation 
whose efficiency will be that of an ordinary triangulation. 

\sa `CGAL::Triangulation_2<Traits,Tds>` 
\sa `CGAL::Delaunay_triangulation_2<Traits,Tds>` 
\sa `TriangulationHierarchyVertexBase_2`
\sa `CGAL::Triangulation_hierarchy_vertex_base_2<Vb>` 

*/
template< typename Tr >
class Triangulation_hierarchy_2 : public Tr {
public:

/// @}

}; /* end Triangulation_hierarchy_2 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgTriangulation2VertexFaceClasses

The class `Triangulation_hierarchy_vertex_base_2` is designed to be used as a vertex base class 
of a triangulation plugged into a 
`Triangulation_hierarchy_2<Tr>`. 

It is 
a model of the concept 
`TriangulationHierarchyVertexBase_2` which refines 
the concept 
`TriangulationVertexBase_2`. 

This class is templated by a parameter `Vb` 
which is to be instantiated by a model of the concept 
`TriangulationVertexBase_2`. 
The class `Triangulation_hierarchy_vertex_base_2<Vb>` inherits 
from the class `Vb`. 
This design allows to use either the default 
vertex base class or a user customized 
vertex base with additional functionalities. 

\cgalModels `TriangulationHierarchyVertexBase_2`

\sa `TriangulationVertexBase_2` 
\sa `TriangulationHierarchyVertexBase_2` 
\sa `CGAL::Triangulation_vertex_base_2<Traits>` 

*/
template< typename Vb >
class Triangulation_hierarchy_vertex_base_2 : public Vb {
public:

/// @}

}; /* end Triangulation_hierarchy_vertex_base_2 */
} /* end namespace CGAL */
