namespace CGAL {

/*!
\ingroup PkgIntervalSkipList

The class `Interval_skip_list` is a dynamic data structure that 
allows to find all members of a set of intervals that overlap a point. 

\cgalHeading{Implementation}

The insertion and deletion of a segment in the interval skip list 
takes expected time \f$ O(\log^2 n)\f$, if the segment endpoints are 
chosen from a continuous distribution. A stabbing query takes expected 
time \f$ O(\log n)\f$, and finding all intervals that contain a point 
takes expected time \f$ O(\log n + k)\f$, where \f$ k\f$ is the number of 
intervals. 

The implementation is based on the code developed by Eric N. Hansen.

*/
template< typename Interval >
class Interval_skip_list {
public:

/// \name Types 
/// @{

/*!
the type of inf and sup of the interval. 
*/ 
typedef Interval::Value Value; 

/*!
An iterator over all intervals. 
*/ 
typedef unspecified_type const_iterator; 

/// @} 

/// \name Creation 
/// @{

/*!
%Default constructor. 
*/ 
Interval_skip_list(); 

/*!
Constructor that inserts the iterator range `[first, last)` in the interval skip list. 
\tparam InputIterator must be an input iterator with value type `Interval`. 
*/ 
template < class InputIterator > 
Interval_skip_list( 
InputIterator first, 
InputIterator last); 

/// @} 

/// \name Operations 
/// @{

/*!
Inserts the iterator range `[first, last)` in the interval skip list, and returns 
the number of inserted intervals. 
\tparam InputIterator must be an input iterator with value type `Interval`. 
*/ 
template < class InputIterator > 
int insert( 
InputIterator first, 
InputIterator last); 

/*!
Inserts interval `i` in the interval skip list. 
*/ 
void insert(const Interval& i); 

/*!
Removes interval `i` from the interval skip list. Returns `true` iff removal was successful. 
*/ 
bool remove(const Interval& i); 

/*!
Returns `true` iff there is an interval that contains `v`. 
*/ 
bool is_contained(const Value& v); 

/*!
Writes the intervals `i` with `i.inf()` \f$ \leq\f$ `v` \f$ \leq\f$ `i.sup` to the 
output iterator `out`. 
\tparam OutputIterator must be an output iterator with value type `Interval`. 
*/ 
template < class OutputIterator > 
OutputIterator find_intervals( 
const Value& v, 
OutputIterator out); 

/*!
Removes all intervals from the interval skip list. 
*/ 
void 
clear(); 

/*!
Returns an iterator over all intervals. 
*/ 
const_iterator begin() const; 

/*!
Returns the past the end iterator. 
*/ 
const_iterator end() const; 

/// @}

}; /* end Interval_skip_list */

/*!
Inserts the interval skip list `isl` into the stream `os`. 
\pre The output operator must be defined for `Interval`. 
\relates Interval_skip_list 
*/ 
template<typenme Interval>
ostream& operator<<(ostream& os, 
const Interval_skip_list<Interval>& isl); 

} /* end namespace CGAL */
