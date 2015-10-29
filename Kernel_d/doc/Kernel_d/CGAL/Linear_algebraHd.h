
namespace CGAL {

/*!
\ingroup PkgKernelDLinAlgClasses

The class `Linear_algebraHd` serves as the default traits class for the LA 
parameter of `CGAL::Homogeneous_d<RT,LA>`. It implements linear 
algebra for Euclidean ring number types `RT`. 

\cgalModels `LinearAlgebraTraits_d`

To make a ring number type `RT` work with this class it has to
provide a division `operator/` with remainder. 

\cgalHeading{Operations}

Fits all operation requirements of the concept. 

*/
template< typename RT >
class Linear_algebraHd {
public:

/// @}

}; /* end Linear_algebraHd */
} /* end namespace CGAL */
