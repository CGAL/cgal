 /// \ingroup PkgSurfaceModeling
 /// \cgalConcept
 ///
 /// @brief Concept describing the set of requirements for quaternion to be used in CGAL::Deform_mesh::rotate

class SurfaceModelingQuaternion
{
public:
  /// Rotate v by quaternion
  /// @tparam Vect a model of SurfaceModelingVect
  /// @param v vector to be rotated
  /// @return rotated vector
  template<class Vect>
  Vect operator*(Vect v);
};