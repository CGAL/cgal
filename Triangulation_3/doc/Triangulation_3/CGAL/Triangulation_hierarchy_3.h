
namespace CGAL {

/*!
\ingroup PkgTriangulation3TriangulationClasses

\deprecated This class is deprecated since \cgal 3.6. Its
functionality has been replaced by the use of the `Fast_location` tag
as the `LocationPolicy` template parameter in
`Delaunay_triangulation_3`.

The class `Triangulation_hierarchy_3` implements a triangulation augmented 
with a data structure which allows fast point location queries. As proved 
in \cite cgal:d-dh-02, this structure has an optimal behavior when it is built 
for Delaunay triangulations. It can however be used for other triangulations. 


\tparam Tr must be instantiated by one of the \cgal triangulation classes. <I>In the current implementation, only 
`Delaunay_triangulation_3` is supported for `Tr`.</I> 
- `Tr::Vertex` has to be a model of the concept `TriangulationHierarchyVertexBase_3`. 
- `Tr::Geom_traits` has to be a model of the concept `DelaunayTriangulationTraits_3`. 

`Triangulation_hierarchy_3` offers exactly the same functionalities as `Tr`. 
Most of them (point location, insertion, removal \f$ \ldots\f$ ) are overloaded to 
improve their efficiency by using the hierarchic structure. 

Note that, since the algorithms that are provided are randomized, the 
running time of constructing a triangulation with a hierarchy may be 
improved when shuffling the data points. 

However, the I/O operations are not overloaded. So, writing a 
hierarchy into a file will lose the hierarchic structure and reading 
it from the file will result in an ordinary triangulation whose 
efficiency will be the same as `Tr`. 

### Implementation ###

The data structure is a hierarchy 
of triangulations. The triangulation at the lowest level is 
the original triangulation where operations and point location are to 
be performed. 
Then at each succeeding level, the data structure 
stores a triangulation of a small random sample of the vertices 
of the triangulation at the preceding level. Point location 
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
data. 

\sa `CGAL::Triangulation_hierarchy_vertex_base_3` 
\sa `CGAL::Delaunay_triangulation_3` 

*/
template< typename Tr >
class Triangulation_hierarchy_3 : public Tr {
public:
}; /* end Triangulation_hierarchy_3 */
} /* end namespace CGAL */
