/*!
\ingroup PkgHDVFConcepts
\cgalConcept

The concept `Ring` describes the requirements for the ring of coefficients used to compute homology in the `CGAL::HDVF` concept. Besides ring operators, it also specifies the functions needed to test invertibility in the ring.
 
\cgalHasModelsBegin
\cgalHasModelsBare{`CGAL::Z`}
\cgalHasModelsBare{`CGAL::Zp`}
\cgalHasModelsEnd

*/

class Ring
{
public:
/// \name Types
/// @{

/*!
 Type of scalars of the ring.
 */
typedef unspecified_type Scalar ;
/// @}
    
    
/// \name Operators
/// @{

/*!
Returns the vector of vertices coordinates.
 */
const std::vector<Point >& get_vertices_coords() const;

/*!
 Returns the coordinates of the cell of index `i` and dimension 0.
 */
Point get_vertex_coords (int i) const;
    
/// @}

};
