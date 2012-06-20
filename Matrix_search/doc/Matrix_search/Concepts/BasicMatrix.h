/// \ingroup PkgMatrixSearchConcepts
///
///
///
/// A class `BasicMatrix` has to provide the following 
/// types and operations in order to be a model for 
/// `BasicMatrix`. 
/// 
/// 
///
///
///
///
/// \hasModel CGAL::Dynamic_matrix<M> 
///
///
/// \sa `SortedMatrixSearchTraits` 
///
class BasicMatrix {
public:

 

 

/// \name Types
/// @{
/*!
\advanced The type of a matrix entry. It has to define
 a copy constructor.
*/
typedef Hidden_type Value;
/// @}

/// \name Operations
/// @{
/*!
\advanced returns the number of columns.
*/
int number_of_columns() const;
/// @}

/// \name Operations
/// @{
/*!
\advanced returns the number of rows.
*/
int number_of_rows() const;
/// @}

/// \name Operations
/// @{
/*!
\advanced returns the entry at position (`row`, `column`). 

\pre \f$ 0 \le\f$ `row` \f$ <\f$ `number_of_rows()`
\pre \f$ 0 \le\f$ `column` \f$ <\f$ `number_of_columns()`
*/
Entry operator()( int row, int column) const;
/// @}

}; /* concept BasicMatrix */

 
 

