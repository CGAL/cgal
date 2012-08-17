
namespace CGAL {

/*!
\ingroup PkgNef2

The class `Filtered_extended_homogeneous` serves as a traits class for the class 
`CGAL::Nef_polyhedron_2<T>`. It uses a polynomial component 
representation based on a ring number type `RT`. 

\models ::ExtendedKernelTraits_2 

Operations 
-------------- 

Fits all operation requirements of the concept. 

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
