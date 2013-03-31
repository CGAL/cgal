 /// \ingroup PkgSurfaceModelingConcepts
 /// \cgalConcept
 ///
 /// @brief Concept describing the set of requirements for quaternion to be used in CGAL::Deform_mesh::rotate

class SurfaceModelingQuaternion
{
public:
/// \name Types 
/// @{
	typedef Hidden_type Vect; ///< a model of SurfaceModelingVect
/// @} 

/// \name Operations 
/// @{
  /// Rotate v by quaternion
  /// @param v vector to be rotated
  /// @return rotated vector
  Vect operator*(Vect v);
/// @}
};