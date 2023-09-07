
/*!
\ingroup PkgCircularKernel3AlgebraicConcepts
\cgalConcept

Concept to represent the roots of a system of three equations of degree 2
in three variables `x`, `y` and `z` that are models of concept
`AlgebraicKernelForSpheres::PolynomialForSpheres_2_3`.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Root_for_spheres_2_3}
\cgalHasModelsEnd

\sa `AlgebraicKernelForSpheres`

*/

class AlgebraicKernelForSpheres::RootForSpheres_2_3 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Test equality of two roots of systems.
*/
bool operator ==(const AlgebraicKernelForSpheres::RootForSpheres_2_3 & r1,
                 const AlgebraicKernelForSpheres::RootForSpheres_2_3 & r2);

/// @}

}; /* end AlgebraicKernelForSpheres::RootForSpheres_2_3 */

