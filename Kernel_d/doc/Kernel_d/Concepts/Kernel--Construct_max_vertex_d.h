
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

\cgalRefines{AdaptableUnaryFunction}

*/
class Kernel_d::Construct_max_vertex_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
returns the lexicographically largest vertex of `ib`.
*/
Kernel_d::Point_d operator()(const Kernel_d::Iso_box_d& ib);


/// @}

};

