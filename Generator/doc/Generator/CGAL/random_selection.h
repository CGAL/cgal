namespace CGAL {

/*!
\ingroup PkgGenerators

chooses a random item from the range \f$ [first,last)\f$ and
writes it to `result`, each item from the range with equal
probability, and repeats this \f$ n\f$ times, thus writing \f$ n\f$ items to
`result`.
A single random number is needed from `rnd` for each item.
Returns the value of `result` after inserting the \f$ n\f$ items.
\pre `Random` is a random number generator type as provided by the STL or by `Random`.

`random_selection` chooses \f$ n\f$ items at random from a random 
access iterator range which is useful to produce degenerate input data 
sets with multiple entries of identical items. 

\sa `CGAL::perturb_points_2` 
*/
template <class RandomAccessIterator, class Size, 
class OutputIterator, class Random>
OutputIterator random_selection( RandomAccessIterator first,
RandomAccessIterator last, 
Size n, OutputIterator result, Random& rnd = default_random);

} /* namespace CGAL */
