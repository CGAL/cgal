
/*!
\ingroup PkgSphericalKernel3Concepts
\cgalconcept


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
/// An object `fo` of this type must provide: 
/// @{

/*! 
Compares the \f$ \theta\f$-coordinates of \f$ p\f$ and \f$ q\f$ in the cylindrical coordinate system relative to the context sphere used by the function `SphericalKernel::compare_theta_3_object`. 
\pre `p` and `q` lie on the context sphere used by the function `SphericalKernel::compare_theta_3_object`, but do not coincide with its poles. 

*/ 
Comparison_result operator() 
(const SphericalKernel::Circular_arc_point_3 &p, 
const SphericalKernel::Circular_arc_point_3 &q ); 

/*! 
Compares the \f$ \theta\f$-coordinates of \f$ p\f$ and of the meridian defined by \f$ m\f$ (see section \ref sectionSKobjects) in the cylindrical coordinate system relative to the context sphere used by the function `SphericalKernel::compare_theta_3_object`. 
\pre `p` lies on the context sphere used by the function `SphericalKernel::compare_theta_3_object`, but does not coincide with its poles. \f$ m \neq(0,0,0)\f$ and the \f$ z\f$-coordinate of \f$ m\f$ is \f$ 0\f$. 

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
Compares the \f$ \theta\f$-coordinates of the meridians defined by \f$ m1\f$ and by \f$ m2\f$ (see section \ref sectionSKobjects) 
in the cylindrical coordinate system relative to the context sphere used by the function `SphericalKernel::compare_theta_3_object`. 
\f$ m1 \neq(0,0,0)\f$, \f$ m2 \neq(0,0,0)\f$ and the \f$ z\f$-coordinate of \f$ m1\f$ and \f$ m2\f$ is \f$ 0\f$. 
*/ 
Comparison_result operator() 
(const SphericalKernel::Vector_3 &m1, 
const SphericalKernel::Vector_3 &m2 ); 

/// @}

}; /* end SphericalKernel::CompareTheta_3 */

