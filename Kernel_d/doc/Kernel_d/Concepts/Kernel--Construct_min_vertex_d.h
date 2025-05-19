
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

\cgalRefines{AdaptableUnaryFunction}

*/
class Kernel_d::Construct_min_vertex_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
returns the lexicographically smallest vertex of `ic`.
*/
Kernel_d::Point_d operator()(const Kernel_d::Iso_box_d& ib);


/// @}

};

