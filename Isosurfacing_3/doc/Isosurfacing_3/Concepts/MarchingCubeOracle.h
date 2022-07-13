/*!
\ingroup PkgIsosurfacing3Concepts
\cgalConcept

The concept `MarchingCubeOracle` describes the set of requirements to be
fulfilled by any class used by the marching cubes algorithms.>`.

\cgalHasModel `CGAL::Marching_cubes_implicit_3`

*/

class MarchingCubeOracle {
public:

/// \name Types
/// @{

/*!
The scalar type.
*/
typedef unspecified_type FT;

/*!
Traits type model of \cgal %Kernel
*/
typedef unspecified_type Traits;


/// @}

};

