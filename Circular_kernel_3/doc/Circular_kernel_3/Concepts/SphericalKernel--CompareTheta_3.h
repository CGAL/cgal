
/*!
\ingroup PkgSphericalKernel3GeometricConcepts
\cgalConcept


\sa `SphericalKernel::CompareThetaZ_3`
\sa `SphericalKernel::CompareX_3`
\sa `SphericalKernel::CompareY_3`
\sa `SphericalKernel::CompareZ_3`
\sa `SphericalKernel::CompareXY_3`
\sa `SphericalKernel::CompareXYZ_3`
\sa `SphericalKernel::Equal_3`

*/

class SphericalKernel::CompareTheta_3 {
public:

/// \name Operations
/// An object of this type must provide: 
/// @{

/*!
Compares the \f$ \theta\f$-coordinates of `p` and `q` in the cylindrical coordinate system relative to the context sphere used by the function `SphericalKernel::compare_theta_3_object`. 
\pre `p` and `q` lie on the context sphere used by the function `SphericalKernel::compare_theta_3_object`, but do not coincide with its poles. 

*/ 
Comparison_result operator() 
(const SphericalKernel::Circular_arc_point_3 &p, 
const SphericalKernel::Circular_arc_point_3 &q ); 

/*!
Compares the \f$ \theta\f$-coordinates of `p` and of the meridian defined by `m` (see Section \ref sectionSKobjects) in the cylindrical coordinate system relative to the context sphere used by the function `SphericalKernel::compare_theta_3_object`. 
\pre `p` lies on the context sphere used by the function `SphericalKernel::compare_theta_3_object`, but does not coincide with its poles. `m` \f$ \neq(0,0,0)\f$ and the \f$ z\f$-coordinate of `m` is \f$ 0\f$. 

*/ 
Comparison_result operator() 
(const SphericalKernel::Circular_arc_point_3 &p, 
const SphericalKernel::Vector_3 &m ); 

/*!
Same as previous, with opposite result. 
*/ 
Comparison_result operator() 
(const SphericalKernel::Vector_3 &m,const SphericalKernel::Circular_arc_point_3 &p); 

/*!
Compares the \f$ \theta\f$-coordinates of the meridians defined by `m1` and by `m2` (see Section \ref sectionSKobjects) 
in the cylindrical coordinate system relative to the context sphere used by the function `SphericalKernel::compare_theta_3_object`. 
`m1` \f$ \neq(0,0,0)\f$, `m2` \f$ \neq(0,0,0)\f$ and the \f$ z\f$-coordinate of `m1` and `m2` is \f$ 0\f$. 
*/ 
Comparison_result operator() 
(const SphericalKernel::Vector_3 &m1, 
const SphericalKernel::Vector_3 &m2 ); 

/// @}

}; /* end SphericalKernel::CompareTheta_3 */

