
namespace CGAL {

/*!
\ingroup PkgNef2

The class `Extended_cartesian` serves as a traits class for the class 
`Nef_polyhedron_2<T>`. It uses a polynomial component 
representation based on a field number type `FT`. 

\cgalModels `ExtendedKernelTraits_2`

### Requirements ###

To make a field number type `FT_model` work with this class, you 
must provide a specialization of the traits class 
`Number_type_traits<FT_model>` 



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
