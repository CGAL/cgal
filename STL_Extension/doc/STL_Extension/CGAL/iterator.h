/// \defgroup STLIterators Iterators and Iterator/Circulator Adaptors
/// \ingroup PkgSTLExtensionRef

namespace CGAL {

/*!
\ingroup STLIterators

The class `Const_oneset_iterator` defines a
`RandomAccessIterator` that always refers to a copy of a
specific object of type `T`.



\cgalModels `RandomAccessIterator`

\sa `CGAL::Emptyset_iterator`
\sa `CGAL::Oneset_iterator`


*/
template< typename T >
class Const_oneset_iterator {
public:

/// \name Creation
/// @{
/*!
creates an iterator that always refers to some copy of
`t`. The copy is constructed by invoking `T`'s copy
constructor and remains constant during the lifetime of the iterator.
*/
Const_oneset_iterator(T&
t);



/// @}



}; /* end Const_oneset_iterator */


/*!
\ingroup STLIterators


The iterator adaptor `Counting_iterator` adds a
counter to the internal iterator of type `Iterator` and defines
equality of two instances in terms of this counter. It can be used
to create finite sequences of possibly infinite sequences of values
from input iterators.



\cgalModels `InputIterator`

\cgalHeading{Requirements}

`Iterator` is a model for
`InputIterator`.

\sa `CGAL::copy_n()`


*/
template< typename Iterator, typename Value >
class Counting_iterator {
public:

/// \name Creation
/// @{
/*!
initializes
the internal counter to `n` and `i` has a singular value.
*/
Counting_iterator( std::size_t n = 0);



/// @}


/// \name Creation
/// @{
/*!

initializes the internal counter to `n` and `i` to `j`.
*/
Counting_iterator( Iterator j, std::size_t n = 0);



/// @}



}; /* end Counting_iterator */


/*!
\ingroup STLIterators



The class `Dispatch_or_drop_output_iterator` defines an
`OutputIterator` that contains a tuple of output iterators, and
dispatches among those based on the type of the value type which is
put in it. Besides defining assignment for all parameters of `V`
and for a tuple of type `V`,  it is also defined for the types `boost::variant<T...>` and
`boost::optional<boost::variant<T...> >`, where `T...`
must be a subset of the parameters of `V`. Should the
`boost::optional` be empty, it will be discarded.

\cgalHeading{Parameters}

\tparam V must be a `std::tuple<...>` of the types of values to be accepted and dispatched.
\tparam O must be a `std::tuple<...>` of the types of corresponding output iterators.

\cgalModels `OutputIterator`

cgalExample{STL_Extension/Dispatch_output_iterator.cpp}

\sa `CGAL::Dispatch_output_iterator<V,O>`
*/
template< typename V, typename O >
class Dispatch_or_drop_output_iterator : public O {
public:

/// \name Types
/// @{
/*!

*/
typedef V Value_type_tuple;



/// @}


/// \name Types
/// @{
/*!

*/
typedef O Iterator_tuple;



/// @}


/// \name Creation
/// @{
/*!
Constructor taking all the output iterators.
*/
Dispatch_or_drop_output_iterator(I...o);

/// @}


/// \name Creation
/// @{
/*!
returns a reference to the tuple of output iterators.
*/
const Iterator_tuple& get_iterator_tuple() const;

/// @}



}; /* end Dispatch_or_drop_output_iterator */

/*!
\returns a `Dispatch_or_drop_output_iterator` constructed from the arguments.
\relates Dispatch_or_drop_output_iterator
*/
template < typename... V, typename... O>
Dispatch_or_drop_output_iterator<tuple<V...>, tuple<O...> >
dispatch_or_drop_output(O... o);



/*!
\ingroup STLIterators



The class `Dispatch_output_iterator` defines an
`OutputIterator` that contains a tuple of output iterators, and
dispatches among those based on the type of the value type which is
put in it. Other types are also accepted, and the object is
discarded in this case. Besides defining assignment for all
parameters of `V` and for a tuple of type `V`, it is also defined for the types
`boost::variant<T...>` and
`boost::optional<boost::variant<T...> >`, where `T...`
can be a list of arbitrary types.

  It also inherits from `O`, which makes it easy to treat like a
  tuple.

\cgalHeading{Parameters}

\tparam V must be a `std::tuple<...>` of the types of values to be accepted and dispatched.
\tparam O must be a `std::tuple<...>` of the types of corresponding output iterators.

\cgalModels `OutputIterator`

\sa `CGAL::Dispatch_or_drop_output_iterator<V,O>`
*/
template< typename V, typename O >
class Dispatch_output_iterator : public O {
public:

/// \name Types
/// @{
/*!

*/
typedef V Value_type_tuple;



/// @}


/// \name Types
/// @{
/*!

*/
typedef O Iterator_tuple;



/// @}


/// \name Creation
/// @{
/*!
Constructor taking all the output iterators.
*/
Dispatch_output_iterator(I...o);



/// @}


/// \name Creation
/// @{
/*!
returns a reference to the tuple of output iterators.
*/
const Iterator_tuple& get_iterator_tuple() const;



/// @}



}; /* end Dispatch_output_iterator */

/*!
\returns a `Dispatch_output_iterator` constructed from the arguments.
\relates Dispatch_output_iterator
*/
template < typename... V, typename... O>
Dispatch_output_iterator<tuple<V...>, tuple<O...> >
dispatch_output(O... o);



/*!
\ingroup STLIterators


The class `Emptyset_iterator` defines an
`OutputIterator` that ignores everything written to it. One can
think of it as being connected to <TT>/dev/null</TT>.



\cgalModels `OutputIterator`

\sa `CGAL::Oneset_iterator`
\sa `CGAL::Const_oneset_iterator`


*/

struct Emptyset_iterator {

/// \name Creation
/// @{
/*!
default constructor.
*/
Emptyset_iterator();



/// @}



}; /* end Emptyset_iterator */


/*!

Constructs `Filter_iterator<Iterator, Predicate>(e, p, c)`.
\relates Filter_iterator
*/
template < class Iterator, class Predicate >
inline Filter_iterator< Iterator, Predicate >
filter_iterator(Iterator e, const Predicate& p, Iterator c = e);

/*!
\ingroup STLIterators



The iterator adaptor `Filter_iterator` acts as a
filter on a given range. Whenever the iterator is in- or
decremented, it ignores all iterators for which the given
`Predicate` is true. The iterator category is the same as for
`Iterator`.

\attention Boost also provides the same functionality via the
`boost::filter_iterator` class. Unfortunately, the semantics
chosen for accepting or rejecting elements based on the predicate's
result are opposite as the semantic chosen here. What is more, the
argument of the predicate is different: the predicate used with
`boost::filter_iterator` must take the value type of the iterator, as
argument, and not the iterator itself.


\tparam Iterator must be a model of `ForwardIterator`
\tparam Predicate must be a functor `Iterator` \f$ \rightarrow\f$  `bool`

*/
template< typename Iterator, typename Predicate >
struct Filter_iterator {

/// \name Creation
/// @{
/*!

*/
Filter_iterator();



/// @}


/// \name Creation
/// @{
/*!
creates an iterator which filters values according to `p`.
Initializes by taking the first valid iterator (according to `p`),
starting at `c`, and stopping at `e` if none is found.
*/
Filter_iterator(Iterator e, Predicate p, Iterator c = e);



/// @}



}; /* end Filter_iterator */


/*!
Constructs `Insert_iterator<Container>(x)`.
\relates Insert_iterator
*/
template < class Container >
Insert_iterator<Container>
inserter(Container &c);

/*!
\ingroup STLIterators



The output iterator `Insert_iterator` is similar
to `std::insert_iterator`, but differs in that it calls the
`insert()` function of the container without the iterator
additional argument.

\cgalModels `OutputIterator`

\tparam Container provides a member function `insert(const Container::const_reference&)`.

*/
template< typename Container >
class Insert_iterator {
public:

/// \name Creation
/// @{
/*!
initializes
the internal container reference to `c`.
*/
Insert_iterator( Container &c );



/// @}



}; /* end Insert_iterator */


/*!
\ingroup STLIterators



The class `Inverse_index` constructs an inverse
index for a given range `[i,j)` of two iterators or circulators of
type `IC`. The first element `I` in the range `[i,j)` has the
index 0. Consecutive elements are numbered incrementally. The
inverse index provides a query for a given iterator or circulator
`k` to retrieve its index number.
\pre The iterator
or circulator must be either of the random access category or the
dereference operator must return stable and distinguishable
addresses for the values, e.g. proxies or non-modifiable iterator
with opaque values will not work.



\cgalHeading{Implementation}

For random access iterators or circulators, it is done in constant
time by subtracting `i`. For other iterator categories, an \stl
`map` is used, which results in a `log(j-i)` query time. The
comparisons are done using the operator `operator<` on pointers.

\sa `CGAL::Random_access_adaptor<IC>`
\sa `CGAL::Random_access_value_adaptor<IC,T>`


*/
template< typename IC >
class Inverse_index {
public:

/// \name Creation
/// @{
/*!
invalid index.
*/
Inverse_index();



/// @}


/// \name Creation
/// @{
/*!
empty inverse
index initialized to start at `i`.
*/
Inverse_index( const IC& i);



/// @}


/// \name Creation
/// @{
/*!
inverse index initialized with range `[i,j)`.
*/
Inverse_index( const IC& i, const IC& j);



/// @}


/// \name Operations
/// @{
/*!
returns inverse index of `k`.
\pre `k` has been stored in the inverse
index.
*/
std::size_t operator[]( const IC& k);



/// @}


/// \name Operations
/// @{
/*!
adds `k` at the end of the indices.
*/
void push_back( const IC& k);



/// @}



}; /* end Inverse_index */


/*!
\ingroup STLIterators



The class `Join_input_iterator_1` joins an iterator and a creator
function object. The result is again an iterator (of the same
iterator category type as the original iterator) that reads an object
from the stream and applies a creator function object to that
object.



\cgalModels `InputIterator`

\sa `CGAL::Creator_1<Arg, Result>`


*/
template< typename Iterator, typename Creator >
class Join_input_iterator_1 {
public:


/// \name Types
/// @{
/*!
*/
typedef Creator::result_type value_type;
/// @}

/// \name Creation
/// @{
/*!
creates a join iterator from the given iterator `i`
and the functor `creator`. Applies `creator` to each item
read from `i`.
*/
Join_input_iterator_1( Iterator i, const Creator&
creator);



/// @}


/// \name Creation
/// @{
/*!
creates a join
iterator from the given iterator `i` and a default constructed
instance of `Creator`. The latter instance is applied to each
item read from `i`.
*/
Join_input_iterator_1( Iterator i);



/// @}



}; /* end Join_input_iterator_1 */


/*!
\ingroup STLIterators

The class `Join_input_iterator_2` joins two iterators. The result is again an iterator (of the same
iterator category type as the original iterator) that reads an object
from the stream and applies a function object to that object.

\cgalModels `InputIterator`


*/
template< typename I1, typename I2, typename Op >
class Join_input_iterator_2 {
public:

  typedef typename Op::result_type value_type;
  typedef typename std::iterator_traits<I1>::difference_type difference_type;
  typedef value_type* pointer;
  typedef value_type& reference;

/// \name Creation
/// @{
/*!
%Default constructor.
 */
Join_input_iterator_2();
/*!
creates a join iterator from the given iterators `i1`, `i2`,
and the functor `op`.
*/
Join_input_iterator_2(I1 i1,I2 i2,const Op& op=Op());

/// @}

/*! returns current position of the first iterator.
 */
 I1 current_iterator1() const { return i1; }

/*! returns current position of the second iterator.
 */
  I2 current_iterator2() const { return i2; }

/*!
  returns `op(current_iterator1(), current_iterator2())`.
 */
  const value_type& operator*() const;



}; /* end Join_input_iterator_2 */

/*!
\ingroup STLIterators

The class `Join_input_iterator_3` joins three iterators. The result is again an iterator (of the same
iterator category type as the original iterator) that reads an object
from the stream and applies a function object to that object.

\cgalModels `InputIterator`


*/
template< typename I1, typename I2,  typename I2, typename Op >
class Join_input_iterator_3 {
public:

  typedef typename Op::result_type value_type;
  typedef typename std::iterator_traits<I1>::difference_type difference_type;
  typedef value_type* pointer;
  typedef value_type& reference;

/// \name Creation
/// @{
/*!
%Default constructor.
 */
Join_input_iterator_3();
/*!
creates a join iterator from the given iterators `i1`, `i2`, `i3`,
and the functor `op`.
*/
  Join_input_iterator_3(I1 i1,I2 i2, I3 i3, const Op& op=Op());

/// @}

/*! returns current position of the first iterator.
 */
 I1 current_iterator1() const { return i1; }

/*! returns current position of the second iterator.
 */
  I2 current_iterator2() const { return i2; }

 /*! returns current position of the second iterator.
 */
  I3 current_iterator3() const { return i3; }

/*!
  returns `op(current_iterator1(), current_iterator2(), current_iterator3())`.
 */
  const value_type& operator*() const;



}; /* end Join_input_iterator_3 */



/*!
\ingroup STLIterators



The adaptor `N_step_adaptor` changes the step width of the
iterator or circulator class `I` to `N`. It is itself an
iterator or circulator respectively. The behavior is undefined if
the adaptor is used on a range `[i,j)` where `j-i` is not a multiple
of `n`.

*/
template< typename I, typename int N >
class N_step_adaptor {
public:

/// \name Creation
/// @{
/*!
down cast.
*/
N_step_adaptor(const I& j);



/// @}



}; /* end N_step_adaptor */


/*!
\ingroup STLIterators



The class `Oneset_iterator` defines an
`BidirectionalIterator` that always refers to one specific
object of type `T`. Internally, `Oneset_iterator` stores a
pointer to the referred object.



\cgalModels `BidirectionalIterator`

\sa `CGAL::Emptyset_iterator`
\sa `CGAL::Const_oneset_iterator`


*/
template< typename T >
class Oneset_iterator {
public:

/// \name Creation
/// @{
/*!
creates
an iterator referring to `t`.
*/
Oneset_iterator(T& t);



/// @}



}; /* end Oneset_iterator */


/*!
\ingroup STLIterators



The class `Random_access_adaptor` provides a random
access for data structures. Either the data structure supports
random access iterators or circulators where this class maps
function calls to the iterator or circulator, or a \stl
`std::vector` is used to provide the random access. The iterator
or circulator of the data structure are of type `IC`.



\sa `CGAL::Inverse_index<IC>`
\sa `CGAL::Random_access_value_adaptor<IC,T>`


*/
template< typename IC >
class Random_access_adaptor {
public:


/// \name Types
/// @{
/*!
size type of the \stl `std::vector`.
*/
typedef unspecified_type size_type;
/// @}

/// \name Creation
/// @{
/*!
invalid index.
*/
Random_access_adaptor();



/// @}


/// \name Creation
/// @{
/*!
empty random
access index initialized to start at `i`.
*/
Random_access_adaptor( const IC& i);



/// @}


/// \name Creation
/// @{
/*!
random access index initialized to the range `[i,j)`.
*/
Random_access_adaptor( const IC& i, const IC& j);



/// @}


/// \name Creation
/// @{
/*!
reserve `r` entries, if a
`std::vector` is used internally.
*/
void
reserve( size_type r);



/// @}


/// \name Operations
/// @{
/*!
returns iterator or circulator to the `n`-th item.
\pre `n <` number of items in the data-structure.
*/
IC operator[]( size_type n);



/// @}


/// \name Operations
/// @{
/*!
adds `k` at the end of the indices.
*/
void push_back( const IC& k);



/// @}



}; /* end Random_access_adaptor */


/*!
\ingroup STLIterators



The class `Random_access_value_adaptor` provides a random
access for data structures. It is derived from
`Random_access_adaptor<IC>`. Instead of returning iterators from
the `operator[]` methods, it returns the dereferenced value of
the iterator. The iterator or circulator of the data structure are
of type `IC`. Their value type is `T`.

\sa `CGAL::Inverse_index<IC>`
\sa `CGAL::Random_access_adaptor<IC>`


*/
template< typename IC, typename T >
class Random_access_value_adaptor {
public:

/// \name Operations
/// Creation and operations see `Random_access_adaptor<IC>`, with
/// the exception of:
/// @{
/*!
returns a reference to the `n`-th item.
\pre `n <` number of items in the data-structure.
*/
T& operator[]( size_type n);



/// @}



}; /* end Random_access_value_adaptor */

} /* end namespace CGAL */
