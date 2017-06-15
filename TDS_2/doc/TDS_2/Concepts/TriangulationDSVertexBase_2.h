
/*!
\ingroup PkgTDS2Concepts
\cgalConcept

The concept `TriangulationDSVertexBase_2` describes the requirements for the 
vertex base class of a `CGAL::Triangulation_data_structure_2<Vb,Fb>`.

Note that if the `CGAL::Triangulation_data_structure_2` 
is plugged into a triangulation class, 
the vertex base class may have additional geometric 
requirements depending on the triangulation class. 

At the base level, 
provides access to one of its incident 
faces through a `Face_handle`. 

Since the `CGAL::Triangulation_data_structure_2` is the class 
which defines the handle 
types, the vertex base class has to be somehow 
parameterized by the triangulation 
data structure. But since the `CGAL::Triangulation_data_structure_2` 
itself is parameterized by the face and vertex 
base classes, there is a cycle in the definition of these classes. 
In order 
to break the cycle, the base classes for faces and vertices 
which are plugged in to instantiate a 
`Triangulation_data_structure_2` 
use  `void` as triangulation 
data structure parameter. Then, 
the `CGAL::Triangulation_data_structure_2` 
uses a <I>rebind</I> mechanism (similar to the one specified in 
`std::allocator`) in order to plug itself 
as parameter in the face and vertex base classes. 
This mechanism requires that the base class provides 
a templated nested class `Rebind_TDS` that 
itself provides 
the subtype `Rebind_TDS::Other` 
which is the <I>rebound</I> version of the base class. 
This <I>rebound</I> base class is the class 
that the `CGAL::Triangulation_data_structure_2` 
actually uses as a base class for the class 
of `CGAL::Triangulation_data_structure_2::Vertex`. 

\cgalRefines `TriangulationDataStructure_2::Vertex` 

\cgalHasModel `CGAL::Triangulation_ds_vertex_base_2<Tds>`
\cgalHasModel `CGAL::Triangulation_vertex_base_2<Traits,Vb>`
\cgalHasModel `CGAL::Regular_triangulation_vertex_base_2<Traits,Vb>`
\cgalHasModel `CGAL::Triangulation_hierarchy_vertex_base_2<Vb>`
\cgalHasModel `CGAL::Triangulation_vertex_base_with_info_2<Info,Traits,Vb>`

\sa `TriangulationVertexBase_2` 
\sa `TriangulationDSFaceBase_2` 
\sa `TriangulationFaceBase_2` 
\sa `TriangulationDataStructure_2::Vertex` 
\sa `CGAL::Triangulation_data_structure_2<Vb,Fb>` 

*/

class TriangulationDSVertexBase_2 {
public:

/// \name Types 
/// The concept `TriangulationDSVertexBase_2` has to provide the following types.
/// @{

/*!
This template class has to define a type `Rebind_TDS<TDS2>::%Other` which is the 
<I>rebound</I> vertex base , where the actual 
`CGAL::Triangulation_data_structure_2` is plugged in. 
This type `Other` will be the actual base 
of the class `CGAL::Triangulation_data_structure_2::Vertex`. 
\note It can be implemented using a nested template class.
*/ 
template <typename TDS2> 
using Rebind_TDS = unspecified_type;

/*!

*/ 
typedef TriangulationDataStructure_2 Triangulation_data_structure; 

/*!

*/ 
typedef TriangulationDataStructure_2::Vertex_handle Vertex_handle; 

/*!

*/ 
typedef TriangulationDataStructure_2::Face_handle Face_handle; 

/// @} 

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
TriangulationDSVertexBase_2(); 

/*!
Constructs a vertex pointing to face `f`. 
*/ 
TriangulationDSVertexBase_2(Face_handle f); 

/// @} 

/// \name Various 
/// These member functions are required by
/// `CGAL::Triangulation_data_structure_2` because it uses
/// `CGAL::Compact_container` to store its faces. See the documentation of
/// `CGAL::Compact_container` for the exact requirements.
/// @{

/*!

*/ 
void * for_compact_container() const; 

/*!

*/ 
void * & for_compact_container(); 

/// @}

}; /* end TriangulationDSVertexBase_2 */

