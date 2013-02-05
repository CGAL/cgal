
/*!
\ingroup PkgPeriodic2Triangulation2Concepts
\cgalConcept

At the base level (see Sections \ref P2Triangulation2secdesign 
and \ref TDS2secdesign), a face stores handles to its four vertices 
and to its four neighbor faces. The vertices and neighbors are 
indexed 0, 1 and 2. Neighbor \f$ i\f$ lies opposite to vertex \f$ i\f$. 

\refines ::TriangulationDSFaceBase_2 

\hasModel CGAL::Periodic_2_triangulation_ds_face_base_2 

\sa `TriangulationDataStructure_2` 
\sa `TriangulationDSFaceBase_2` 
\sa `Periodic_2TriangulationDSVertexBase_2` 

*/

class Periodic_2TriangulationDSFaceBase_2 {
public:

/// \name Access Functions 
/// @{

/*! 
Returns the offset of vertex `i`. 
\pre \f$ i \in\{0, 1, 2\}\f$. 
*/ 
int offset(int i) const; 

/*! 
Returns true if the offset of vertex `i` is zero for \f$ i \in\{0, 1, 2\}\f$. 
*/ 
bool has_zero_offsets() const; 

/// @} 

/// \name Setting 
/// @{

/*! 
Sets the vertex offsets according to `off0` to `off2`. 
*/ 
void set_offsets(int off0, int off1, int off2); 

/// @}

}; /* end Periodic_2TriangulationDSFaceBase_2 */

