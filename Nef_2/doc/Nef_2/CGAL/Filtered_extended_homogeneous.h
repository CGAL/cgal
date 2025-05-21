
namespace CGAL {

/*!
\ingroup PkgNef2Ref

The class `Filtered_extended_homogeneous` serves as a traits class for the class
`Nef_polyhedron_2<T>`. It uses a polynomial component
representation based on a ring number type `RT`.

\cgalModels{ExtendedKernelTraits_2}

\tparam RT  must be a model of `RingNumberType`.

\sa `CGAL::Extended_cartesian<FT>`
\sa `CGAL::Extended_homogeneous<RT>`

*/
template< typename RT >
class Filtered_extended_homogeneous {
public:

/// \name Creation
/// @{

/*!
default constructor.
*/
Filtered_extended_homogeneous();

/// @}

}; /* end Filtered_extended_homogeneous */
} /* end namespace CGAL */
