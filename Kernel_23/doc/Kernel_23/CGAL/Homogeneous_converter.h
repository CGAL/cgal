
namespace CGAL {

/*!
\ingroup kernel_conversion

`Homogeneous_converter` converts objects from the kernel traits `K1` to
the kernel traits `K2`. Those traits must be of the form
`Homogeneous<RT1>` and `Homogeneous<RT2>` (or the equivalent with
`Simple_homogeneous`). It then provides the following operators to
convert objects from `K1` to `K2`.

The third template parameter `RTConverter` is a function object that must
provide `K2::RT operator()(const K1::RT &n);` that
converts `n` to an `K2::RT` that has the same value.

The default value of this parameter is `CGAL::NT_converter<K1::RT, K2::RT>`,
which uses the conversion operator from
`K1::RT` to `K2::RT`.

Similarly, the fourth template parameter must provide
`K2::FT operator()(const K1::FT &n);` that
converts `n` to an `K2::FT` that has the same value. Its
default value is `CGAL::NT_converter<K1::FT, K2::FT>`.

\sa `CGAL::Homogeneous<RingNumberType>`
\sa `CGAL::Simple_homogeneous<RingNumberType>`

*/
template< typename K1, typename K2, typename RTConverter, typename FTConverter >
class Homogeneous_converter {
public:

/// \name Creation
/// @{

/*!
%Default constructor.
*/
Homogeneous_converter<>();

/// @}

/// \name Operations
/// Similar operators as the one shown are defined for the other
/// kernel traits geometric types `Point_3`, `Vector_2`...
/// @{

/*!
returns a `K2::Point_2` which coordinates are those of `p`,
converted by `RTConverter`.
*/
K2::Point_2 operator()(const K1::Point_2&p);

/// @}

}; /* end Homogeneous_converter */
} /* end namespace CGAL */
