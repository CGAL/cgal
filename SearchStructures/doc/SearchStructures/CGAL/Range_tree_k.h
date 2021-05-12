
namespace CGAL {

/*!
\ingroup PkgSearchStructuresSearchStructure

An object of the class `Range_tree_k` is a \f$ k\f$-dimensional range tree
that can store k-dimensional keys of type `Key`.
The class allows to perform
window queries on the keys. The class `Range_tree_k` is parameterized with
a range tree traits class `Traits` that defines, among other things,
the type of the `Key`.

\cgal provides traits class implementations that allow to use
the range tree with point classes from the \cgal kernel as keys.
These classes are `CGAL::Range_segment_tree_traits_set_2<R>`,
`CGAL::Range_segment_tree_traits_set_3<R>`,
`CGAL::Range_tree_traits_map_2<R>` and
`CGAL::Range_tree_traits_map_3<R>`. The concept
RangeSegmentTreeTraits_d defines the requirements that range tree traits
classes must fulfill. This allows the advanced user to develop further
range tree traits classes.

\cgalHeading{Example}

The following example program uses the predefined `Range_tree_2` data structure together with the predefined traits
class `Range_tree_map_traits_2` which has two template
arguments specifying the
type of the point data in each dimension
(`CGAL::Cartesian<double>`) and the value type of the
2-dimensional point data (`char`). Therefore the `Range_tree_2` is defined on 2-dimensional point data
(`CGAL::Point_2<Cartesian<double> >`) each of which is
associated with a character.
Then, a few data items are created and put into a list. After
that the tree is constructed according to that list, a window
query is performed, and the query elements are given out.

\code
#include <CGAL/Cartesian.h>
#include <CGAL/Range_segment_tree_traits.h>
#include <CGAL/Range_tree_k.h>

typedef CGAL::Cartesian<double> K;
typedef CGAL::Range_tree_map_traits_2<K, char> Traits;
typedef CGAL::Range_tree_2<Traits> Range_tree_2_type;

int main()
{
  typedef Traits::Key Key;
  typedef Traits::Interval Interval;

  std::vector<Key> InputList, OutputList;
  InputList.push_back(Key(K::Point_2(8,5.1), 'a'));
  InputList.push_back(Key(K::Point_2(1,1.1), 'b'));
  InputList.push_back(Key(K::Point_2(3,2.1), 'c'));

  Range_tree_2_type Range_tree_2(InputList.begin(),InputList.end());
  Interval win(Interval(K::Point_2(4,8.1), K::Point_2(5,8.2)));
  std::cout << "\n Window Query:\n ";
  Range_tree_2.window_query(win, std::back_inserter(OutputList));
  std::vector<Key>::iterator current=OutputList.begin();
  while(current!=OutputList.end()){
    std::cout << (*current).first.x() << "," << (*current).first.y()
              << ":" << (*current++).second << std::endl;
  }
}
\endcode
*/
template< typename Traits >
class Range_tree_k {
public:

/// \name Types
/// @{

/*!
the type of the range tree traits class.
*/
typedef unspecified_type Traits;

/*!

*/
typedef Traits::Key Key;

/*!

*/
typedef Traits::Interval Interval;

/// @}

/// \name Creation
/// @{

/*!
Introduces an empty range tree `R`.
*/
Range_tree_k ();

/*!
Introduces a range tree `R` and initializes it with the data
in the range `[first, last)`.
\pre `value_type(first) == Traits::Key`.
*/
template < class ForwardIterator >
Range_tree_k (ForwardIterator first,
ForwardIterator last);

/// @}

/// \name Operations
/// @{

/*!
Introduces a range tree `R` and initializes it with the data
in the range `[first, last)`. This function can only be applied
once on an empty range tree.
\pre `value_type(first) == Traits::Key`.
*/
template < class ForwardIterator >
void
make_tree(ForwardIterator first,
ForwardIterator last);

/*!
writes all data that are in the interval `window` to the container
where `out` points to, and returns an output iterator that points
to the last location the function wrote to.
\pre `value_type(out) == Traits::Key`.
*/
template < class OutputIterator >
OutputIterator
window_query(Interval window,
OutputIterator out);

/// @}

}; /* end Range_tree_k */
} /* end namespace CGAL */
