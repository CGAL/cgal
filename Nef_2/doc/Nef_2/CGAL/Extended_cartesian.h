
namespace CGAL {

/*!
\ingroup PkgNef2Ref

The class `Extended_cartesian` serves as a traits class for the class
`Nef_polyhedron_2<T>`. It uses a polynomial component
representation based on a field number type `FT`.

\cgalModels{ExtendedKernelTraits_2}

\tparam FT  must be a model of `FieldNumberType`.



\sa `CGAL::Extended_homogeneous<RT>`
\sa `CGAL::Filtered_extended_homogeneous<RT>`

*/
template< typename FT >
class Extended_cartesian {
public:

/// \name Creation
/// @{

/*!
default constructor.
*/
Extended_cartesian();

/// @}

}; /* end Extended_cartesian */
} /* end namespace CGAL */
