namespace CGAL {
/// \ingroup PkgMatrixSearch
///
///
///
///
/// The class `Sorted_matrix_search_traits_adaptor` can be used as an 
/// adaptor to create sorted matrix search traits classes for 
/// arbitrary feasibility test and matrix classes `F` resp. 
/// `M`. 
/// 
///
///
///
/// \models ::SortedMatrixSearchTraits 
///
///
/// Requirements 
/// -------------- 
/// 
/// <OL> 
/// <LI>`M` is a model for `BasicMatrix` <I>and</I> 
/// <LI>`F` defines a copy constructor and a monotone `bool 
/// operator()( const Value&)`. 
/// </OL> 
/// 
/// 
/// 
///
///
template< F,M >
class Sorted_matrix_search_traits_adaptor {
public:

 

 

/// \name Creation
/// @{
/*!
\advanced initializes `t` to use `m` for feasibility testing.
*/
Sorted_matrix_search_traits_adaptor<F,M>( const F&
 m);
/// @}

/// \name Types
/// @{
/*!
\advanced typedef to `M`.
*/
typedef Hidden_type Matrix;
/// @}

/// \name Types
/// @{
/*!
\advanced typedef to `Matrix::Value`.
*/
typedef Hidden_type Value;
/// @}

/// \name Types
/// @{
/*!
\advanced typedef to
 `std::less<Value>`.
*/
typedef Hidden_type Compare_strictly;
/// @}

/// \name Types
/// @{
/*!
\advanced typedef to
 `std::less_equal<Value>`.
*/
typedef Hidden_type Compare_non_strictly;
/// @}

/// \name Operations
/// @{
/*!
\advanced returns the `Compare_strictly` object to be used for the search.
*/
Compare_strictly compare_strictly()
 const;
/// @}

/// \name Operations
/// @{
/*!
\advanced returns the `Compare_non_strictly` object to be used for the search.
*/
Compare_non_strictly compare_non_strictly()
 const;
/// @}

/// \name Operations
/// @{
/*!
\advanced uses the feasibility test given during creation.
*/
bool is_feasible(const Value& a);
/// @}

 

 
}; /* class Sorted_matrix_search_traits_adaptor */
} /* namespace CGAL */

 
 

