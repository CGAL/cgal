
namespace CGAL {

/*!
\ingroup PkgRangeSegmentTreesDSearchStructure

An object of the class `Segment_tree_k` is a \f$ k\f$-dimensional segment tree 
that can store k-dimensional intervals of type `Interval`. 
The class allows to perform 
window queries, enclosing queries, and inverse range queries on the keys. The class `Segment_tree_k` is parameterized with 
a segment tree traits class `Traits` that defines, among other things, 
the type of the `Interval`. 
In order to perform an inverse range query, a range query of 
\f$ \epsilon\f$ width has to be performed. We prefered not to offer an 
extra function for this sort of query, since the inverse range 
query is a special case of the range query. Furthermore, offering 
an inverse range query in the segment tree class implies offering this 
function also in the range tree class and having an extra item in 
the traits class that accesses the inverse range query point. 

\cgal provides traits class implementations that allow to use 
the segment tree with point classes from the \cgal kernel as keys. 
These classes are `CGAL::Range_segment_tree_traits_set_2<R>`, 
`CGAL::Range_segment_tree_traits_set_3<R>`, 
`CGAL::Segment_tree_traits_map_2<R>` and 
`CGAL::Segment_tree_traits_map_3<R>`. The concept 
RangeSegmentTreeTraits_d defines the requirements that segment tree traits 
classes must fulfill. This allows the advanced user to develop further 
segment tree traits classes. 

\cgalHeading{Example}

This example illustrates the use of the predefined segment tree 
on 3-dimensional interval data (with no value associated). After 
the definition of the traits type and tree type, some intervals 
are constructed and the tree is build according to the 
intervals. Then, a window query is performed and the query 
elements are given out. 

\code
#include <CGAL/Cartesian.h> 
#include <CGAL/Segment_tree_k.h> 
#include <CGAL/Range_segment_tree_traits.h> 

typedef CGAL::Cartesian<int> K; 
typedef CGAL::Range_segment_tree_set_traits_3<K> Traits; 
typedef CGAL::Segment_tree_3<Traits> Segment_tree_3_type; 

int main() 
{ 
  typedef Traits::Interval Interval; 
  typedef Traits::Key Key; 
  std::list<Interval> InputList, OutputList; 

  InputList.push_back(Interval(Key(1,5,7), Key(2,7,9))); 
  InputList.push_back(Interval(Key(2,7,6), Key(3,8,9))); 
  InputList.push_back(Interval(Key(6,9,5), Key(9,13,8))); 
  InputList.push_back(Interval(Key(1,3,4), Key(3,9,8))); 

  Segment_tree_3_type Segment_tree_3(InputList.begin(),InputList.end()); 

  Interval a(Key(3,6,5), Key(7,12,8)); 
  Segment_tree_3.window_query(a,std::back_inserter(OutputList)); 
  std::list<Interval>::iterator j = OutputList1.begin(); 
  std::cout << "\n window_query (3,6,5),(7,12,8) \n"; 
  while(j!=OutputList.end()){ 
    std::cout << (*j).first.x() << "," << (*j).first.y() << ","; 
    std::cout << (*j).first.z() <<", " << (*j).second.x() << ","; 
    std::cout << (*j).second.y() << "," << (*j).second.z() << std::endl; 
    j++; 
  } 
} 
\endcode
*/
template< typename Traits >
class Segment_tree_k {
public:

/// \name Types 
/// @{

/*!
the type of the segment tree traits class. 
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
Introduces an empty segment tree `S`. 
*/ 
Segment_tree_k (); 

/*!
Introduces a segment tree `S` and initializes it with the data 
in the range `[first, last)`. 
\pre `value_type(first) == Traits::Interval`. 
*/ 
template < class ForwardIterator > 
Segment_tree_k (ForwardIterator first, 
ForwardIterator last); 

/// @} 

/// \name Operations 
/// @{

/*!
Introduces a segment tree `S` and initializes it with the data 
in the range `[first, last)`. This function can only be applied 
once on an empty segment tree. 
\pre `value_type(first) == Traits::Interval`. 
*/ 
template < class ForwardIterator > 
void 
make_tree(ForwardIterator first, 
ForwardIterator last); 

/*!
writes all intervals that have non empty intersection with interval `window` to the container 
where `out` points to, and returns an output iterator that points 
to the last location the function wrote to. 
\pre `value_type(out) == Traits::Interval`. 
*/ 
template < class OutputIterator > 
OutputIterator 
window_query(Interval window, 
OutputIterator out); 

/*!
writes all intervals that enclose in the interval `window` to the container 
where `out` points to, and returns an output iterator that points 
to the last location the function wrote to. 
\pre `value_type(out) == Traits::Interval`. 
*/ 
template < class OutputIterator > 
OutputIterator 
enclosing_query(Interval window, 
OutputIterator out); 

/// @}

}; /* end Segment_tree_k */
} /* end namespace CGAL */
