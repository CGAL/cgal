
CONVERROR Additional namespace TriangulationDataStructure:: required
/*!
\ingroup PkgTriangulationsConcepts
\cgalconcept

The concept `Vertex` represents the vertex class of a triangulation 
data structure. It must define 
the types and operations listed in this section. Some of these 
requirements are of geometric nature, they are <I>optional</I> 
when using the triangulation data structure class alone. They become 
compulsory when the triangulation data structure is used as a layer 
for the geometric triangulation class. 

Creation 
-------------- 

In order to obtain new vertices or destruct unused vertices, the user must 
call the `new_vertex()` or `delete_vertex()` method of the 
triangulation data structure. 

\sa `TriangulationDataStructure::FullCell` 
CONVERRORSeeAlso: `TriangulationDataStructure`. 

*/

class Vertex {
public:

/// \name Types 
CONVERROR Check if this needs to be spread\n/// The class `Vertex` defines types that are the same as some of the types defined by the triangulation data structure class `TriangulationDataStructure`.
/// @{

/*! 
<I>Optional for the triangulation data 
structure alone.</I> 
*/ 
typedef Hidden_type Point; 

/*! 

*/ 
typedef TriangulationDataStructure Triangulation_data_structure; 

/*! 

*/ 
typedef TriangulationDataStructure::Vertex_handle Vertex_handle; 

/*! 

*/ 
typedef TriangulationDataStructure::Full_cell_handle Full_cell_handle; 

/// @} 

/// \name Access Functions 
/// @{

/*! 
Returns a full cell of the triangulation having `v` as vertex. 
*/ 
Full_cell_handle full_cell() const; 

/*! 
Returns the point stored in the vertex. 
<I>Optional for the triangulation data structure alone.</I> 
*/ 
Point point() const; 

/// @} 

/// \name Setting 
CONVERROR Check if this needs to be spread\n/// CONVERROR DEBUG
/// @{

/*! 
Sets the incident cell to `c`. 
*/ 
void set_full_cell(Full_cell_handle c); 

/*! 
Sets the point to `p`. <I>Optional for the 
triangulation data structure alone.</I> 
*/ 
void set_point(const Point & p); 

/// @} 

/// \name Checking 
/// @{

/*! 
Checks the validity of the vertex. Must check that its incident cell 
has this vertex. The validity of the base vertex is also checked. 

When `verbose` is set to `true`, messages are printed to give 
a precise indication on the kind of invalidity encountered. 
*/ 
bool is_valid(bool verbose = false) const; 

/// @}

}; /* end Vertex */

