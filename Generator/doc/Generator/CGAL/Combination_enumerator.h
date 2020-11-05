namespace CGAL {

/// \ingroup PkgGeneratorsRef

/*!
  The class `Combination_enumerator` is used to enumerate all fixed-size combinations
  (subsets) of a <i>source range</i> of elements. For example, it can
  enumerate all the combinations of 2 elements from the source range `[3,7)`
  (`7` excluded) which gives the enumeration `{3,4}, {3,5}, {3,6}, {4,5},
  {4,6}, {5,6}`. The source range consists of elements of type
  `CombinationElement` and is specified by its first element and the element
  just beyond its last one.

  \tparam CombinationElement should be a model of the concept `CombinationElement`.

  Each combination is uniquely represented as an increasing sequence of elements.
  Thus, the combinations can be lexicographically ordered. They are enumerated in
  that order, so that we can talk about the first or last combination.
*/
template <class CombinationElement>
class Combination_enumerator{
public:
/// \name Creation
/// @{

/*!
  This constructor initializes the object to
enumerate the combinations of `k` elements from the source range
`[first, beyond)`. The current combination is set to the first combination
of the enumeration.
\pre  `1 <= k <= beyond - first`
*/
Combination_enumerator(int k, const CombinationElement & first, const CombinationElement & beyond);

/*!
The copy constructor.
*/
Combination_enumerator(const Combination_enumerator & combi);

/// @}

/// \name Access to the Current Combination
/// @{

/*!
Returns the `i`-th element of the current combination.
\pre `0 <= i < number_of_elements()`
*/
const CombinationElement & operator[](int i);

/// @}


/// \name Access to the Enumeration
/// @{

/*!
Returns the size of the enumerated combinations (the parameter `k` from the class' constructor).
*/
int number_of_elements();

/*!
Returns the smallest element of the source range. (the parameter `first` of the
constructor of the class).
*/
const CombinationElement & min_element();

/*!
Returns the successor to the largest element of the source range (the parameter `beyond` of
the constructor of the class).
*/
const CombinationElement & beyond_element();


/*!
Returns `true` if and only if all combinations have been enumerated.
*/
bool finished();


/// @}


/// \name Operations
/// @{

/*!
Resets the enumerator. The current combination is set to the first one of the enumeration.
*/
void reset();

/*!
Moves *this to the next combination.
*/
void operator++();


/*!
Post-incrementation. Same as the pre-incrementation above, but returns the original value of `*this`.
*/
Combination_enumerator operator++(int);
/// @}

}; /*end of class Combination_enumerator*/
} /* end namespace CGAL */
