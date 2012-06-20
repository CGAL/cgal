namespace CGAL {
/// \ingroup PkgMatrixSearch
///
///
///
///
/// The class `Dynamic_matrix` is an adaptor for an 
/// arbitrary matrix class `M` to provide the dynamic operations 
/// needed for monotone matrix search. 
/// 
///
///
///
/// Requirements 
/// -------------- 
/// `M` is a model for `BasicMatrix`. 
/// 
/// 
///
///
///
/// \models ::BasicMatrix 
///
///
/// \sa `BasicMatrix` 
///
///
/// Implementation 
/// -------------- 
/// All operations take constant time except for 
/// `extract_all_even_rows` which needs time linear in the number 
/// of rows. 
/// 
/// 
/// 
/// 
///
///
template< M >
class Dynamic_matrix {
public:

 

 

/// \name Creation
/// @{
/*!
\advanced initializes `d` to `m`. `m` is <I>not</I> copied, we only store a reference.
*/
Dynamic_matrix( const M& m);
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

/// \name Operations
/// @{
/*!
\advanced replace column `old` with column number `new`. 
\pre \f$ 0 \le\f$ `old`, `new` \f$ <\f$ `number_of_columns()`
*/
void replace_column( int old, int new);
/// @}

/// \name Operations
/// @{
/*!
\advanced returns a new Matrix consisting of all rows of `d` with even index, (i.e. first row is row \f$ 0\f$ of `d`, second row is row \f$ 2\f$ of `d` etc.). \pre `number_of_rows()` \f$ > 0\f$.
*/
Matrix* extract_all_even_rows() const;
/// @}

/// \name Operations
/// @{
/*!
\advanced deletes the rightmost columns, such that `d` becomes quadratic. 

\pre `number_of_columns()` \f$ \ge\f$ `number_of_rows()`
\post `number_of_rows()` \f$ ==\f$ `number_of_columns()`
*/
void shrink_to_quadratic_size();
/// @}

}; /* class Dynamic_matrix */
} /* namespace CGAL */

 
 

