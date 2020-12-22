
/*!
\ingroup PkgCircularKernel3AlgebraicConcepts
\cgalConcept

*/

class AlgebraicKernelForSpheres::Solve {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Copies in the output iterator the common roots of `p1`, `p2`,
and `p3`, with their multiplicity, as objects of type
`std::pair< AlgebraicKernelForSpheres::Root_for_spheres_2_3, int>`.

Here, `Type1`, `Type2`, and `Type3` can all be either
`AlgebraicKernelForSpheres::Polynomial_1_3` or
`AlgebraicKernelForSpheres::Polynomial_for_spheres_2_3`.

\pre The set of solutions of the system is 0-dimensional.
*/
template < class OutputIterator >
OutputIterator
operator()(const Type1 &p1,
const Type2 &p2,
const Type3 &p3,
OutputIterator res);

/// @}

}; /* end AlgebraicKernelForSpheres::Solve */

