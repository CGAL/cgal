
namespace CGAL {

/*!
\ingroup PkgNef2

The class `Extended_homogeneous` serves as a traits class for the class 
`Nef_polyhedron_2<T>`. It uses a polynomial component 
representation based on a Euclidean ring number type `RT`. 

\cgalModels `ExtendedKernelTraits_2`

### Requirements ###

To make an Euclidean ring number type 
`RT_model` work with this class the number type must support 
a gcd computation in namespace `CGAL::NTS`. \cgal provides 
a function template for this, which will be used by default when 
your number type is not one of the built-in number types, one of 
the number types distrubuted with \cgal or one of the \leda 
number types. 


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
