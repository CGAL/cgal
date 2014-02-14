/// \ingroup PkgSurfaceModelingConcepts
/// \cgalConcept
///
/// Concept describing the set of requirements of a simple point type.
///
class RawPoint_3
{
public:
/// \name Creation
/// @{
  RawPoint_3();
  RawPoint_3(double x, double y, double z);
/// @}

/// \name Coordinates Accessors
/// @{
  double x() const;
  double y() const;
  double z() const;
  /// `i<=0 && i<3`
  double& operator[](int i);
  /// `i<=0 && i<3`
  double  operator[](int i) const;
/// @}
};
