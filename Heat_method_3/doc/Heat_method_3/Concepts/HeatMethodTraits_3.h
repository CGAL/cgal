/*!
\ingroup PkgHeatMethodConcepts

\cgalConcept

The concept `HeatMethodTraits_3` describes the types,
predicates, and constructions required by the traits class parameter of
`CGAL::Heat_method_3::Heat_method_3`.

\cgalHasModel `CGAL::Heat_method::Heat_method_Eigen_traits_3`.


*/

class HeatMethodTraits_3
{
public:

/// \name Types
/// @{

  /// Mesh type, is required for Heat_method_3
  typedef unspecified_type Triangle_mesh;

  // compute...
  compute();

  // solve ....
  solve();

/// @}


};
