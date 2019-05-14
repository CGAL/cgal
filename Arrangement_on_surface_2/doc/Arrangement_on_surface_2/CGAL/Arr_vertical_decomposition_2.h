namespace CGAL {

/*!
\ingroup PkgArrangement2Funcs

Produces the symbolic vertical decomposition of a 
given arrangement, performing a batched vertical ray-shooting query from 
all arrangement vertices, such that every vertex is associated with a pair 
of objects, one corresponds to the arrangement feature that lies below it, 
and the other corresponds to the feature that lies above it. 
The output of this function can be readily used for inserting vertical walls 
and physically decomposing the arrangement into pseudo-trapezoids. To do 
this, it is convenient to process the vertices in an ascending 
\f$ xy\f$-lexicographic order. The visible objects are therefore returned through 
an output iterator, which pairs each finite arrangement vertex with the two 
features it "sees", such that the vertices are given in ascending 
\f$ xy\f$-lexicographic order. 

Produces the symbolic vertical decomposition of the `arr` arrangement. 
More precisely, it performs a batched vertical ray-shooting query from all 
arrangement vertices, such that every vertex is associated with a pair of 
objects, one corresponding to the arrangement feature that lies below it, 
while the other corresponds to the feature that lies above it. 
The query results are returned through the output iterator, which pairs 
each finite arrangement vertex with a pair of `Object`s, the first 
represents the feature below the vertex, and the second represents the 
feature that lies above it. Each `Object` may be one of the following: 
<UL> 
<LI>`Halfedge_const_handle`, if the vertex is located above (or 
below) an edge. The given halfedge is always directed from right to left. 
In case there is no concrete edge below (or above) the vertex, and 
the arrangement is unbounded, then the object returned is a 
<I>fictitious</I> halfedge. 
<LI>`Face_const_handle`, in case there is no edge below (or above) 
the vertex, and the arrangement is bounded. 
<LI>`Vertex_const_handle`, in case the vertex is located vertically 
above (or below) another arrangement vertex. 
<LI>An empty object, in case the vertex is the top end-vertex of 
a vertical edge, we define there is no feature below it. Similarly, if 
it is the bottom end-vertex of a vertical edge, we define that there 
is no feature above it. 
</UL> 
The function returns a past-the-end iterator for its output sequence. 

\cgalHeading{Requirements}

`OutputIterator::value_type` must be 
`pair<Arrangement_2::Vertex_const_handle, pair<Object, Object> >`. 

*/
template<typename Traits, typename Dcel,
typename OutputIterator>
OutputIterator decompose (const Arrangement_2<Traits,Dcel>& arr,
OutputIterator oi);

} /* namespace CGAL */

