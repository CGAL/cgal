namespace CGAL {

/*!
\ingroup PkgGenerators
\brief chooses `n` items at random from a random 
access iterator range which is useful to produce degenerate input data 
sets with multiple entries of identical items. 

The function chooses a random item from the range `[first,last)` and
writes it to `result`, each item from the range with equal
probability, and repeats this \f$ n\f$ times, thus writing `n` items to
`result`.
A single random number is needed from `rnd` for each item.
Returns the value of `result` after inserting the `n` items.
\pre `Random` is a random number generator type as provided by the STL or by `Random`.



\sa `CGAL::perturb_points_2()` 
*/
template <class RandomAccessIterator, class Size, 
class OutputIterator, class Random>
OutputIterator random_selection( RandomAccessIterator first,
RandomAccessIterator last, 
Size n, OutputIterator result, Random& rnd = get_default_random());

} /* namespace CGAL */
