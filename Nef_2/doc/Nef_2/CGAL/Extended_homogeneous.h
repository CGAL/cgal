
namespace CGAL {

/*!
\ingroup PkgNef2Ref

The class `Extended_homogeneous` serves as a traits class for the class
`Nef_polyhedron_2<T>`. It uses a polynomial component
representation based on a Euclidean ring number type `RT`.

\cgalModels `ExtendedKernelTraits_2`


\tparam RT  must be a model of `RingNumberType`.


\sa `CGAL::Extended_cartesian<FT>`
\sa `CGAL::Filtered_extended_homogeneous<RT>`

*/
template< typename RT >
class Extended_homogeneous {
public:

/// \name Creation
/// @{

/*!
default constructor.
*/
Extended_homogeneous();

/// @}

}; /* end Extended_homogeneous */
} /* end namespace CGAL */
