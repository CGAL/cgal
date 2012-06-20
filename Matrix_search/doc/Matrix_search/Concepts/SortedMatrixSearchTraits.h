/// \ingroup PkgMatrixSearchConcepts
///
///
///
///
/// The concept `SortedMatrixSearchTraits` defines types and operations 
/// needed to compute the smallest entry in a set of sorted matrices 
/// that fulfills a certain feasibility criterion using the function 
/// `sorted_matrix_search`. 
/// 
/// 
///
///
///
///
/// \hasModel CGAL::Sorted_matrix_search_traits_adaptor<F,M> 
///
///
/// \sa `BasicMatrix` 
///
class SortedMatrixSearchTraits {
public:

 

 

/// \name Types
/// @{
/*!
\advanced The class used for representing matrices.
 It has to be a model for `BasicMatrix`.
*/
typedef Hidden_type Matrix;
/// @}

/// \name Types
/// @{
/*!
\advanced The class used for representing the matrix elements.
*/
typedef Matrix::Value Value;
/// @}

/// \name Types
/// @{
/*!
\advanced An adaptable binary function
 class: `Value` \f$ \times\f$ `Value` \f$ \rightarrow\f$ `bool`
 defining a non-reflexive total order on `Value`. This
 determines the direction of the search.
*/
typedef Hidden_type Compare_strictly;
/// @}

/// \name Types
/// @{
/*!
\advanced An adaptable binary function
 class: `Value` \f$ \times\f$ `Value` \f$ \rightarrow\f$ `bool`
 defining the reflexive total order on `Value` corresponding
 to `Compare_strictly`.
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
\advanced The predicate to determine whether an element `a` is feasible. It has to be monotone in the sense that `compare( a, b)` and `is_feasible( a)` imply `is_feasible( b)`.
*/
bool is_feasible( const Value& a);
/// @}

}; /* concept SortedMatrixSearchTraits */

 
 

