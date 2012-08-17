
namespace CGAL {

/*!
\ingroup PkgNef2

The class `Extended_cartesian` serves as a traits class for the class 
`CGAL::Nef_polyhedron_2<T>`. It uses a polynomial component 
representation based on a field number type `FT`. 

\models ::ExtendedKernelTraits_2 

Requirements 
-------------- 

To make a field number type `FT_model` work with this class, you 
must provide a traits class for this number type: 
`CGAL::Number_type_traits<FT_model>` (See the support library manual.) 

Operations 
-------------- 

Fits all operation requirements of the concept. 

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
