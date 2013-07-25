/*!
\ingroup PkgIntervalSkipListConcepts
\cgalConcept

The concept `Interval` describes the requirements for the 
template argument `Interval` of a `Interval_skip_list<Interval>`. 

The concept does not specify, whether the interval is open or 
closed. It is up to the implementer of a model for this concept 
to define that. 

\cgalHasModel `CGAL::Interval_skip_list_interval<Value>`
\cgalHasModel `CGAL::Level_interval` 

\sa `Interval_skip_list` 

*/

class Interval {
public:

/// \name Creation 
/// @{

/*!
Default constructor. 
*/ 
Interval(); 

/// @} 

/// \name Types 
/// @{

/*!
The type of the lower and upper bound of the interval. 
*/ 
typedef unspecified_type Value; 

/// @} 

/// \name Access Functions 
/// @{

/*!
Returns the lower bound. 
*/ 
Value inf() const; 

/*!
Returns the upper bound. 
*/ 
Value sup() const; 

/*!
Returns `true`, iff the interval contains `v`. 
*/ 
bool contains(const Value& v) const; 

/*!
Returns `true`, iff the interval contains `(i,s)`. 
*/ 
bool contains_interval(const Value& i, const Value& s) const; 

/*!
Equality test. 
*/ 
bool operator==(const Interval& I) const; 

/*!
Unequality test. 
*/ 
bool operator!=(const Interval& I) const; 

/// @}

}; /* end Interval */
