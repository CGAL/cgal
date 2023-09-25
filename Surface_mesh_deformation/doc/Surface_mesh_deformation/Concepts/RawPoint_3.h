/// \ingroup PkgSurfaceMeshDeformationConcepts
/// \cgalConcept
/// \cgalRefines{DefaultConstructible,Assignable}
///
/// Concept describing the set of requirements of a simple point type.
///
class RawPoint_3
{
public:
/// \name Creation
/// @{
  /// constructor from %Cartesian coordinates
  RawPoint_3(double x, double y, double z);
/// @}

/// \name Coordinates Accessors
/// @{
  /// `i<=0 && i<3`
  double  operator[](int i) const;
/// @}
};
