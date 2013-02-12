
/*!
\ingroup PkgPeriodic2Triangulation2Concepts
\cgalConcept

At the base level (see Section \ref
Section_2D_Triangulations_Software_Design), a face stores handles to
its four vertices and to its four neighbor faces. The vertices and
neighbors are indexed 0, 1 and 2. Neighbor \f$ i\f$ lies opposite to
vertex \f$ i\f$.

\cgalRefines ::TriangulationDSFaceBase_2 

\cgalHasModel CGAL::Periodic_2_triangulation_ds_face_base_2 

\sa `TriangulationDataStructure_2` 
\sa `TriangulationDSFaceBase_2` 
\sa `Periodic_2TriangulationDSVertexBase_2` 

*/

class Periodic_2TriangulationDSFaceBase_2 {
public:

/// \name Access Functions 
/// @{

// TODO(NGHK): Check
/*! 
Returns the offset of vertex `i`. 
\pre \f$ i \in\{0, 1, 2\}\f$. 
*/ 
int offset(int i) const; 

// TODO(NGHK): Check
/*! 
Returns true if the offset of vertex `i` is zero for \f$ i \in\{0, 1, 2\}\f$. 
*/ 
bool has_zero_offsets() const; 

/// @} 

/// \name Setting 
/// @{

// TODO(NGHK): Check
/*! 
Sets the vertex offsets according to `off0` to `off2`. 
*/ 
void set_offsets(int off0, int off1, int off2); 

/// @}

}; /* end Periodic_2TriangulationDSFaceBase_2 */

