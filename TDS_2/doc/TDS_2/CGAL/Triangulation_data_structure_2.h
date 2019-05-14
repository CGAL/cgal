
namespace CGAL {

/*!
\ingroup PkgTDS2

The class `Triangulation_data_structure_2` is a model 
for the `TriangulationDataStructure_2` concept. 
It can be used to represent an orientable 2D triangulation 
embedded in a space of any dimension. 

The vertices and faces are stored in two nested containers, which are 
implemented using `Compact_container`. The class may offer some 
flexibility for the choice of container in the future, in the form of 
additional template parameters. 

\cgalModels `TriangulationDataStructure_2`

\cgalHeading{Modifiers}

In addition to the modifiers required by the 
`TriangulationDataStructure_2` concept, the `Triangulation_data_structure_2` class 
supports also the modifiers below. Note also that the modifiers below 
guarantee the combinatorial validity of the resulting data structure. 

\cgalHeading{Illustrations}

\anchor figtdssplitjoin 
\image html join_split.png "The join and split operations."
\image latex join_split.png "The join and split operations."

\anchor figtdsirdeg2 
\image html tds-insert_degree_2.png "Insertion and removal of degree 2 vertices. "
\image latex tds-insert_degree_2.png "Insertion and removal of degree 2 vertices. "
*/
template< typename Vb, typename Fb >
class Triangulation_data_structure_2 {
public:


/// \name Types

/// @{

  typedef Triangulation_data_structure_2<Vb,Fb>  Tds;

/// The vertex type.
  typedef  typename Vb::template Rebind_TDS<Tds>::Other  Vertex;

/// The face type.
  typedef  typename Fb::template Rebind_TDS<Tds>::Other  Face;

/// @}

/// \name  
/// \cgalAdvancedBegin
/// In addition to the interface documented in the concept,
/// the class offers the following types.
/// \cgalAdvancedEnd
/// @{

/*!
Vertex container type. 
*/ 
typedef Compact_container<Vertex> Vertex_range; 

/*!
Face container type. 
*/ 
typedef Compact_container<Face> Face_range; 

/// @} 

/// \name Operations 
/// \cgalAdvancedBegin
/// In addition to the interface documented in the concept,
/// the class offers the following functions.
/// \cgalAdvancedEnd
/// @{

/*!
Returns a reference to the container of faces. 
*/ 
Face_range & faces() const; 

/*!
Returns a reference to the container of faces. 
*/ 
Face_range & faces(); 

/*!
Returns a reference to the container of vertices. 
*/ 
Vertex_range & vertices() const; 

/*!
Returns a reference to the container of vertices. 
*/ 
Vertex_range & vertices(); 

/// @} 

/// \name Modifiers 
/// @{

/*!
Joins 
the vertices that are endpoints of the edge `(f,i)`, and returns 
a vertex handle to common vertex (see 
Fig.\ \ref figtdssplitjoin). 
\pre `f` must be different from `Face_handle()` and `i` must be `0`, `1` or `2`. 
*/ 
Vertex_handle join_vertices(Face_handle f, int i); 

/*!
Joins 
the vertices that are endpoints of the edge `e`, and returns 
a vertex handle to common vertex. 
*/ 
Vertex_handle join_vertices(Edge e); 

/*!
Joins 
the vertices that are endpoints of the edge `*eit`, and returns 
a vertex handle to common vertex. 
*/ 
Vertex_handle join_vertices(Edge_iterator eit); 

/*!
Joins 
the vertices that are endpoints of the edge `*ec`, and returns 
a vertex handle to common vertex. 
*/ 
Vertex_handle join_vertices(Edges_circulator ec); 

/*!
Splits the vertex `v` into two vertices `v1` and 
`v2`. The common faces `f` and `g` of `v1` and 
`v2` are created after (in the counter-clockwise sense) the 
faces `f1` and `f2`. The 4-tuple `(v1,v2,f,g)` is 
returned (see Fig. \ref figtdssplitjoin). 
\pre `dimension()` must be equal to `2`, `f1` and `f2` must be different from `Face_handle()` and `v` must be a vertex of both `f1` and `f2`. 
*/ 
boost::tuples::tuple<Vertex_handle, Vertex_handle, Face_handle, 
Face_handle> 
split_vertex(Vertex_handle v, Face_handle f1, Face_handle 
f2); 

/*!
Inserts 
a degree two vertex and two faces adjacent to it that have two common 
edges. The edge defined by the face handle `f` and the integer 
`i` is duplicated. It returns a handle to the vertex created 
(see Fig. \ref figtdsirdeg2). 
*/ 
Vertex_handle insert_degree_2(Face_handle f, int i); 

/*!
Removes a degree 2 
vertex and the two faces adjacent to it. The two edges of the star of 
`v` that are not incident to it are collapsed 
(see Fig. \ref figtdsirdeg2). 
\pre The degree of `v` must be equal to 2. 
*/ 
void remove_degree_2(Vertex_handle v); 

/// @}

}; /* end Triangulation_data_structure_2 */
} /* end namespace CGAL */
