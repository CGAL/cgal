
/*!
\ingroup PkgPeriodic2Triangulation2Concepts
\cgalConcept

A refinement of the concept `TriangulationDSVertexBase_2` 
which adds an API for offset.

At the base level of 2D-triangulations (see Section \ref
Section_2D_Triangulations_Software_Design), a vertex provides access
to one of its incident faces through a handle.

The storage of the offset is only needed when a triangulation is copied.

\cgalRefines `TriangulationDSVertexBase_2` 

 

\cgalHasModel CGAL::Periodic_2_triangulation_ds_vertex_base_2 

\sa `TriangulationDataStructure_2` 
\sa `TriangulationDSVertexBase_2` 
\sa `Periodic_2TriangulationDSFaceBase_2` 

*/

class Periodic_2TriangulationDSVertexBase_2 {
public:

/// \name Types 
/// @{

/*! 
A model of the concept 
`Periodic_2Offset_2` 
*/ 
typedef Hidden_type Periodic_2_offset_2; 

/// @} 

/// \name Access Functions 
/// @{

/*! 
Returns the offset stored in the vertex. 
*/ 
Periodic_2_offset_2 offset() const; 

/*! 
Returns `true` if the offset has been set, `false` otherwise. 
*/ 
bool get_offset_flag() const; 

/*! 
Sets the offset and sets the offset flag to `true`. 
*/ 
void set_offset(Periodic_2_offset_2 o); 

/*! 
Sets the offset flag to `false` and clears the offset. 
*/ 
void clear_offset(); 

/// @}

}; /* end Periodic_2TriangulationDSVertexBase_2 */

