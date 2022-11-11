namespace CGAL {

/*!
\ingroup AdvancedClasses

The class `Kd_tree_rectangle` implements `d`-dimensional iso-rectangles and related operations,
e.g., methods to compute bounding boxes of point sets.

*/
template< typename FT , typename Dimension>
class Kd_tree_rectangle {
public:

/// \name Types
/// @{

/*!
Dimension type. Either `CGAL::Dimension_tag`
or `CGAL::Dynamic_dimension_tag`.
*/
typedef Dimension Dimension;

/*!
Number type.
*/
typedef FT FT;

/// @}

/// \name Creation
/// @{

/*!
Constructs a `d`-dimensional rectangle `r` with lower bound and upper bound set to zero
in each dimension.
*/
Kd_tree_rectangle(int d);

/*!
Constructs the bounding box of the points in the range `[begin,end)`, where the value
type of `PointIter` can be used by operators of functors `Construct_cartesian_const_iterator_d`
to define iterators with value type `FT`.
*/
template <class Construct_cartesian_const_iterator_d,class PointIter>
Kd_tree_rectangle(int d, PointIter begin, PointIter end,const Construct_cartesian_const_iterator_d& construct_it);

/// @}

/// \name Operations
/// @{

/*!
Returns the lower bound of the rectangle in dimension `i`.
*/
FT min_coord(int i) const;

/*!
Returns the upper bound of the rectangle in dimension `i`.
*/
FT max_coord(int i) const;

/*!
Sets upper bound in dimension `i` to `x`.
*/
void set_upper_bound(int i, const FT& x);

/*!
Sets lower bound in dimension `i` to `x`.
*/
void set_lower_bound(int i, const FT& x);

/*!
Returns the maximal span of the rectangle.
*/
FT max_span() const;

/*!
Returns the smallest coordinate for which the rectangle has its maximal span.
*/
FT max_span_coord() const;

/*!
Returns the dimension of the rectangle.
*/
int dimension() const;

/*!
Splits rectangle in dimension `d` at coordinate-value `value`
by modifying itself to lower half and by modifying `r` to upper half.
*/
void split(Kd_tree_rectangle<FT,Dimension>& r, int d, FT value);

/// @}

}; /* end Kd_tree_rectangle */

/*!
Inserts rectangle `r` in the output stream `s` and returns `s`.
\relates Kd_tree_rectangle
*/
template<class FT>
std::ostream& operator<<(std::ostream& s, Kd_tree_rectangle<FT>& r);

} /* end namespace CGAL */
