namespace CGAL {

/*!
\ingroup PkgSpatialSorting

The class `Multiscale_sort` represents a sorting algorithm adaptor. 
Given a range of \f$ n\f$ points: 
<OL> 
<LI>it applies `Sort` on the last \f$ (1 -\f$`ratio`\f$ )\times n\f$ points, 
<LI>it recurses on the first `ratio`\f$ \times n\f$ points, 
stopping when there are less than `threshold` points. 
</OL> 

*/
template< typename Sort >
class Multiscale_sort {
public:

/// \name Creation 
/// @{

/*! 
constructs an instance with `traits` as traits class instance. 
*/ 
Multiscale_sort (const Sort &sort = Sort(), std::ptrdiff_t threshold = 1, double ratio = 0.5); 

/// @} 

/// \name Operations 
/// @{

/*! 
sorts the range [`begin`, `end`). 
\require `Sort::operator()(RandomAccessIterator begin, RandomAccessIterator end)` is defined. 
*/ 
template <class RandomAccessIterator> void operator() (RandomAccessIterator begin, RandomAccessIterator end) const; 

/// @}

}; /* end Multiscale_sort */
} /* end namespace CGAL */
