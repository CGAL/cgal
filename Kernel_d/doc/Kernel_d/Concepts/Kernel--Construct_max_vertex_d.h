
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

\cgalRefines `AdaptableFunctor` (with one argument) 

*/
class Kernel_d::Construct_max_vertex_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
returns the lexicographically largest vertex of `ic`. 
*/ 
Kernel_d::Point_d operator()(const Kernel_d::Iso_cuboid_d& ic); 


/// @}

};

