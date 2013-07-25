
/*!
\ingroup PkgHandlesAndCirculatorsConcepts
\cgalConcept

A constant iterator range. Refer to the `Range` concept for more details. 

\cgalRefines Boost's Range concept 

\cgalHasModel STL containers 
\cgalHasModel <A HREF="http://www.boost.org/libs/range/doc/html/range/reference/utilities/iterator_range.html">`boost::iterator_range`</A> 

\sa `Range` 

*/
class ConstRange {
public:

/// \name Types 
/// @{

/*!
The constant iterator type. 
*/ 
typedef unspecified_type const_iterator; 

/*!
An unsigned integral type that can represent the 
size of a range. 
*/ 
typedef unspecified_type size_type; 

/// @} 

/// \name Member functions 
/// @{

/*!
returns the const iterator pointing to the first element. 
*/ 
const_iterator begin() const; 

/*!
returns the past-the-end const iterator. 
*/ 
const_iterator end() const; 

/*!
returns the size of the range. 
*/ 
size_type size() const; 

/*!
returns whether the range is empty. 
*/ 
bool empty() const; 

/// @}

}; /* end ConstRange */

