
namespace CGAL {

/*!
\ingroup kernel_conversion

`Cartesian_converter` converts objects from the kernel traits `K1` to 
the kernel traits `K2` using `NTConverter` to do the conversion. 
Those traits must be of the form 
`Cartesian<FT1>` and `Cartesian<FT2>` (or the equivalent with 
`Simple_cartesian`). It then provides the following operators to convert 
objects from `K1` to `K2`. 

The third template parameter `NTConverter` is a function object that must 
provide `K2::FT operator()(K1::FT n)` that converts `n` to an 
`K2::FT` which has the same value. 

The default value of this parameter is `CGAL::NT_converter<K1::FT, K2::FT>`. 

\cgalHeading{Example}

In the following example, we compute exactly 
the intersection point between a line and a triangle, 
and we then create a double approximation of this point. 

\cgalExample{Kernel_23/cartesian_converter.cpp} 

\sa `CGAL::Cartesian<FieldNumberType>` 
\sa `CGAL::Simple_cartesian<FieldNumberType>` 

*/
template< typename K1, typename K2, typename NTConverter >
class Cartesian_converter {
public:

/// \name Creation 
/// @{

/*!
%Default constructor. 
*/ 
Cartesian_converter<>(); 

/// @} 

/// \name Operations 
/// Similar operators are defined for the other kernel traits types `Point_3`, `Vector_2`...
/// @{

/*!
returns a `K2::Point_2` which coordinates are those of `p`, 
converted by `NTConverter`. 
*/ 
K2::Point_2 operator()(const K1::Point_2&p); 

/// @}

}; /* end Cartesian_converter */
} /* end namespace CGAL */
