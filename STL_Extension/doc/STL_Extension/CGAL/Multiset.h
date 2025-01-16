
namespace CGAL {

/*!
\ingroup PkgSTLExtensionRef


\anchor classcgal_multimap


An instance `s` of the parametrized data type `Multiset` is a
multi-set of elements of type `Type`, represented as a red-black tree
(see [\cgalCite{clrs-ia-01} Chapter 13 for an excellent introduction to red-black
trees).
The main difference between `Multiset` and the \stl `std::multiset` is that
the latter uses a less-than functor with a Boolean return type, while our
`Multiset` class is parameterized by a comparison functor `Compare` that
returns the three-valued `Comparison_result` (namely it returns either
`SMALLER`, `EQUAL`, or `LARGER`). It is thus possible to maintain
the underlying red-black tree with less invocations of the comparison functor.
This leads to a speedup of about 5% even if we maintain a set of integers.
When each comparison of two elements of type `Type` is an expensive
operation (for example, when they are geometric entities represented using
exact arithmetic), the usage of a three-valued comparison functor can lead to
considerable decrease in the running times.

Moreover, `Multiset` allows the insertion of an element into the set given
its <I>exact</I> position, and not just using an insertion hint, as done by
`std::multiset`. This can further reduce the running times, as additional
comparison operations can be avoided.

In addition, the `Multiset` guarantees that the order of elements sent to the
comparison functor is fixed. For example, if we insert a new element `x`
into the set (or erase an element from the set), then we always invoke
`Compare() (x, y)` (and never `Compare() (y, x)`), where `y` is an
element already stored in the set. This behavior, not supported by
`std::multiset`, is sometimes crucial for designing more efficient
comparison predicates.

`Multiset` also allows for look-up of keys whose type may differ from
`Type`, as long as users supply a comparison functor `CompareKey`,
where `CompareKey() (key, y)` returns the three-valued
`Comparison_result` (`key` is the look-up key and `y` is an
element of type `Type`). Indeed, it is very convenient to look-up
equivalent objects in the set given just by their key. We note however that
it is also possible to use a key of type `Type` and to employ the default
`Compare` functor for the look-up, as done when using the
`std::multiset` class.

\warning Finally, `Multiset` introduces the `catenate()` and `split()`
functions. The first function operates on `s` and accepts a second
set `s2`, such that the maximum element in `s` is not greater than
the minimal element in `s2`, and concatenates `s2` to `s`. The
second function splits `s` into two sets, one containing all the
elements that are less than a given key, and the other contains all
elements greater than (or equal to) this key.



\tparam Type the type of the stored elements.
\tparam Compare  the comparison-functor type. This type should provide
the following operator for comparing two `Type` elements, namely:
<br>
`Comparison_result operator() (const Type& t1, const Type& t2) const;`
<br>
The `CGAL::Compare<Type>` functor is used by default. In this case,
`Type` must support an equality operator (`operator==`) and a
less-than operator (`operator<`).
\tparam Allocator the allocator type. `CGAL_ALLOCATOR` is used by default.


\cgalHeading{Implementation}

`Multiset` uses a proprietary implementation of a red-black tree
data-structure. The red-black tree invariants guarantee that the height of a
tree containing \f$ n\f$ elements is \cgalBigO{\log{n}} (more precisely, it is bounded by
\f$ 2 \log_{2}{n}\f$). As a consequence, all methods that accept an element and need
to locate it in the tree (namely `insert(x)`, `erase(x)`,
`find(x)`, `count(x)`, `lower_bound(x)` , `upper_bound(x)`,
`find_lower(x)` and `equal_range(x)`) take \cgalBigO{\log{n}} time and
perform \cgalBigO{\log{n}} comparison operations.

On the other hand, the set operations that accept a position iterator (namely
`insert_before(pos, x)`, `insert_after(pos, x)` and `erase(pos)`)
are much more efficient as they can be performed at a <I>constant</I> amortized
cost (see \cgalCite{gs-dfbt-78} and \cgalCite{t-dsna-83} for more details).
More important, these set operations require <I>no</I> comparison operations.
Therefore, it is highly recommended to maintain the set via iterators
to the stored elements, whenever possible. The function `insert(pos, x)`
is safer to use, but it takes amortized \cgalBigO{\min\{d,\log{n}\}} time, where \f$ d\f$
is the distance between the given position and the true position of `x`.
In addition, it always performs at least two comparison operations.

The `catenate()` and `split()` functions are also very efficient, and
can be performed in \cgalBigO{\log{n}} time, where \f$ n\f$ is the total number of
elements in the sets, and without performing any comparison operations
(see \cgalCite{t-dsna-83} for the details).
Note however that the size of two sets resulting from a split operation is
initially unknown, as it is impossible to compute it in less than linear time.
Thus, the first invocation of `size()` on such a set takes linear time,
and <I>not</I> constant time.

The design is derived from the \stl `multiset` class-template (see,
e.g, \cgalCite{cgal:ms-strg-96}), where the main differences between the two
classes are highlighted in the class definition above.

*/
template< typename Type, typename Compare, typename Allocator >
class Multiset {
public:



/*!
\name Types
In compliance with \stl, the types `value_type` and `key_type`
(both equivalent to `Type`), `reference` and `const_reference`
(reference to a value-type), `key_compare` and `value_compare`
(both equivalent to `Compare`), `size_type` and `difference_type`
are defined as well.
*/
/// @{
/*!

*/
typedef unspecified_type iterator;

/*!
bi-directional iterators for the elements stored in the set.
*/
typedef unspecified_type const_iterator;


/*!

*/
typedef unspecified_type reverse_iterator;

/*!
reverse bi-directional iterators for the elements stored in the set.
*/
typedef unspecified_type const_reverse_iterator;

/// @}

/// \name Creation
/// @{
/*!

creates an an empty set `s` that uses a default comparison
functor.
*/
Multiset<Type,Compare,Allocator>();



/// @}


/// \name Creation
/// @{
/*!

creates an an empty set `s` that uses the given comparison
functor `comp`.
*/
Multiset<Type,Compare,Allocator>(
const Compare& comp);



/// @}


/// \name Creation
/// @{
/*!

creates a set `s` containing all elements in the range
`[first, last)`, that uses the comparison
functor `comp`.
*/
template <class InputIterator>
Multiset<Type,Compare,Allocator>(
InputIterator first, InputIterator last,
const Compare& comp = Compare());



/// @}


/// \name Creation
/// @{
/*!

copy constructor.
*/
Multiset<Type,Compare,Allocator>(
const Multiset<Type,Compare,Allocator>& other);



/// @}


/// \name Creation
/// @{
/*!

assignment operator.
*/
const Multiset<Type,Compare,Allocator>& operator= (
const Multiset<Type,Compare,Allocator>& other);



/// @}


/// \name Creation
/// @{
/*!

swaps the contents of `s` with those of the other set.
*/
void swap (Multiset<Type,Compare,Allocator>& other);



/// @}


/// \name Access Member Functions
/// @{
/*!

the comparison functor used.
*/
Compare key_comp() const;



/// @}


/// \name Access Member Functions
/// @{
/*!

the comparison functor used (same as above).
Both functions have a non-const version that return a reference
to the comparison functor.
*/
Compare value_comp() const;



/// @}


/// \name Access Member Functions
/// @{
/*!

returns `true` if the set is empty, `false` otherwise.
*/
bool empty ();



/// @}


/// \name Access Member Functions
/// @{
/*!

returns the number of elements stored in the set.
*/
size_t size ();



/// @}


/// \name Access Member Functions
/// @{
/*!

returns the maximal number of elements the set can store
(same as `size()`).
*/
size_t max_size ();



/// @}


/// \name Access Member Functions
/// @{
/*!

returns an iterator pointing to the first element stored in the set
(a `const` version is also available).
*/
iterator begin();



/// @}


/// \name Access Member Functions
/// @{
/*!

returns an iterator pointing beyond the last element stored in the set
(a `const` version is also available).
*/
iterator end();



/// @}


/// \name Access Member Functions
/// @{
/*!

returns a reverse iterator pointing beyond the last element stored in the
set (a `const` version is also available).
*/
reverse_iterator rbegin();



/// @}


/// \name Access Member Functions
/// @{
/*!

returns a reverse iterator pointing to the first element stored in the set
(a `const` version is also available).
*/
reverse_iterator rend();



/// @}


/// \name Comparison Operations
/// @{
/*!

returns `true` if the sequences of elements in the two sets are
pairwise equal (using the comparison functor).
*/
bool operator== (
const Multiset<Type,Compare,Allocator>& other) const;



/// @}


/// \name Comparison Operations
/// @{
/*!

returns `true` if the element sequence in `s` is
lexicographically smaller than the element sequence of `other`.
*/
bool operator< (
const Multiset<Type,Compare,Allocator>& other) const;



/// @}


/// \name Insertion Methods
/// @{
/*!

inserts the element `x` into the set and returns an iterator pointing
to the newly inserted element.
*/
iterator insert (const Type& x);



/// @}


/// \name Insertion Methods
/// @{
/*!

inserts all elements in the range `[first, last)` into
the set.
*/
template <class InputIterator>
void insert (InputIterator first, InputIterator last);



/// @}


/// \name Insertion Methods
/// @{
/*!

inserts the element `x` with a given iterator used as a hint for the
position of the new element. It Returns an iterator pointing to the
newly inserted element.
*/
iterator insert (iterator position, const Type& x);



/// @}


/// \name Insertion Methods
/// @{
/*!

inserts the element `x` as the predecessor of the element at the given
position.
\pre The operation does not violate the set order - that is,
`x` is not greater than the element pointed by
`position` and not less than its current predecessor.
*/
iterator insert_before (iterator position, const Type& x);



/// @}


/// \name Insertion Methods
/// @{
/*!

inserts the element `x` as the successor of the element at the given
position.
\pre The operation does not violate the set order - that is,
`x` is not less than the element pointed by
`position` and not greater than its current successor.
*/
iterator insert_after (iterator position, const Type& x);



/// @}


/// \name Removal Methods
/// @{
/*!

erases all elements equivalent to `x` from the set and returns the
number of erased elements.
*/
size_t erase (const Type& x);



/// @}


/// \name Removal Methods
/// @{
/*!

erases the element pointed by `position`.
*/
void erase (iterator position);



/// @}


/// \name Removal Methods
/// @{
/*!

clears the set (erases all stored elements).
*/
void clear ();



/// @}


/*!
\name Look-up Methods
All methods listed in this section can also accept a `Type` element
as a look-up key. In this case, it is not necessary to supply a `CompareKey`
functor, as the `Compare` functor will be used by default.
*/
/// @{

/*!

searches for the an element equivalent to `key` in the set. If the
set contains objects equivalent to `key`, it returns an iterator
pointing to the first one. Otherwise, `end()` is returned (a
`const` version is also available).
*/
template <class Key, class CompareKey>
iterator find (const Key& key, const CompareKey& comp_key);

/*!

returns the number of elements equivalent to `key` in the set.
*/
template <class Key, class CompareKey>
size_t count (const Key& key, const CompareKey& comp_key) const;

/*!

returns an iterator pointing to the first element in the set that is not
less than `key`. If all set elements are less than `key`,
`end()` is returned (a `const` version is also available).
*/
template <class Key, class CompareKey>
iterator lower_bound (const Key& key, const CompareKey& comp_key);

/*!

returns an iterator pointing to the first element in the set that is
greater than `key`. If no set element is greater than `key`,
`end()` is returned (a `const` version is also available).
*/
template <class Key, class CompareKey>
iterator upper_bound (const Key& key, const CompareKey& comp_key);

/*!

returns the range of set elements equivalent to the given key, namely
`(lower_bound(key), upper_bound(key))` (a `const` version is
also available).
*/
template <class Key, class CompareKey>
std::pair<iterator,iterator>
equal_range (const Key& key, const CompareKey& comp_key);


/*!

returns a pair comprised of `lower_bound(key)` and a Boolean flag
indicating whether this iterator points to an element equivalent to
the given key (a `const` version is also available).
*/
template <class Key, class CompareKey>
std::pair<iterator,bool>
find_lower (const Key& key, const CompareKey& comp_key);


/// @}


/// \name Special Operations
/// @{
/*!

replaces the element stored at the given position with `x`.
\pre The operation does not violate the set order - that is,
`x` is not less that `position`'s predecessor and
not greater than its successor.
*/
void replace (iterator position, const Type& x);



/// @}


/// \name Special Operations
/// @{
/*!

swaps places between the two elements given by `pos1` and `pos2`.
\pre The operation does not violate the set order - that is,
`pos1` and `pos2` store equivalent elements.
*/
void swap (iterator pos1, iterator pos2);



/// @}


/// \name Special Operations
/// @{
/*!

concatenates all elements in `s_prime` into `s` and clears
`s_prime`.
All iterators to `s` and to `s_prime` remain valid.
\pre The maximal element in `s` is not greater than the minimal
element in `s_prime`.
*/
void catenate (Self& s_prime);



/// @}


/// \name Special Operations
/// @{
/*!

splits `s` such that it contains all elements that are less than
the given `key` and such that `s_prime` contains all other elements.
\pre `s_prime` is initially empty.
*/
template <class Key, class CompareKey>
void split (Key key, CompareKey comp_key, Self& s_prime);



/// @}


/// \name Special Operations
/// @{
/*!

splits `s` such that it contains all set elements in the range
`[begin, position)` and such that `s_prime` contains all elements
in the range `[position, end())`.
\pre `s_prime` is initially empty.
*/
void split (iterator position, Self& s_prime);



/// @}



}; /* end Multiset */
} /* end namespace CGAL */
