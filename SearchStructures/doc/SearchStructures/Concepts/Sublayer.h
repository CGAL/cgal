
/*!
\ingroup PkgRangeSegmentTreesDConcepts
\cgalConcept

Defines the requirements that a 
`Sublayer_type` of class `CGAL::Range_tree_d` 
or `CGAL::Segment_tree_d` has to fulfill. 

First of all, the class has to be derived from the abstract base 
class `CGAL::Tree_base` and therefore 
has to provide methods 
`make_tree`, `window_query`, 
`enclosing_query` and 
`is_inside` 
with the same parameter types as the instantiated class 
`Range_tree_d` or `Segment_tree_d`, 
respectively. 
Furthermore a method `bool is_anchor()` has to be provided. If the `Sublayer_type` class 
builds a recursion anchor for class 
`Segment_tree_d`, this function is expected to 
return `true`, `false` otherwise. 

Such a recursion anchor class is provided by the class class. 
`CGAL::Tree_anchor<Data, Window>`. 

*/

class Sublayer {
public:

/// @}

}; /* end Sublayer */

