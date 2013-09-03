namespace CGAL {

/*!
 * \addtogroup PkgPointSet2RangeSearch
 *
 * There are six versions of the function template `range_search()` that 
 * perform range searches on Delaunay triangulations. The first performs 
 * circular range searches, the second triangular range searches and the 
 * third performs iso-rectangular range searches. The other three range search 
 * function templates perform enhanced variants of the three aforementioned 
 * operations. 
 *
 * They get a user-defined object that has to control the range
 * search operation.  This way one can for instance stop the search,
 * when `n` points were found.
 */

/*!
\ingroup PkgPointSet2RangeSearch

computes handles to all vertices contained in the closure of disk `C`.
The computed vertex handles will be placed as a sequence of objects in a container of value type
of `res`
which points to the first object in the sequence. The function
returns an output iterator pointing to the position beyond the end
of the sequence.
`delau`$ is the \cgal Delaunay triangulation on which we perform the range search operation.

\cgalHeading{Requirements}

<UL> 
<LI>`Dt` is a \cgal Delaunay triangulation and contains the following subset of types from the concept `PointSetTraits` and from 
the Delaunay triangulation data type: 
<UL> 
<LI>`Dt::Geom_traits` 
<LI>`Dt::Vertex_handle` 
<LI>`Dt::Vertex` 
<LI>`Dt::Vertex_circulator` 
<LI>`Dt::Vertex_iterator` 
<LI>`Dt::Point` 
<LI>`Dt::Geom_traits::Bounded_side_2` 
<LI>`Dt::Geom_traits::Construct_center_2` 
</UL> 
<LI>the template parameter `Circle` corresponds to `Dt::Geom_traits::Cricle_2` 
</UL> 

*/
template<class Dt, class Circle, class OutputIterator>
OutputIterator range_search(Dt& delau, const Circle& C, OutputIterator res);

/*!
\ingroup PkgPointSet2RangeSearch

computes handles to all vertices contained in the closure of the triangle `(a,b,c)`.

\pre `a`, `b`, and `c` must not be collinear. 
The computed vertex handles will be placed as a sequence of objects in a container of value type
of `res`
which points to the first object in the sequence. The function
returns an output iterator pointing to the position beyond the end
of the sequence. 
`delau` is the \cgal Delaunay triangulation on which we perform the range search operation.

\cgalHeading{Requirements}

`Dt` is a \cgal Delaunay triangulation and contains the following subset of types from the concept `PointSetTraits` and from 
the Delaunay triangulation data type: 
<UL> 
<LI>`Dt::Geom_traits` 
<LI>`Dt::Vertex_handle` 
<LI>`Dt::Vertex` 
<LI>`Dt::Vertex_circulator` 
<LI>`Dt::Vertex_iterator` 
<LI>`Dt::Point` 
<LI>`Dt::Geom_traits::Circle_2` 
<LI>`Dt::Geom_traits::Bounded_side_2` 
<LI>`Dt::Geom_traits::Construct_center_2` 
<LI>`Dt::Geom_traits::Orientation_2` 
<LI>`Dt::Geom_traits::Construct_circle_2` 
</UL> 

*/
template<class Dt, class OutputIterator>
OutputIterator range_search(Dt& delau, const Dt::Point& a, const Dt::Point& b, 
const Dt::Point& c,OutputIterator res);

/*!
\ingroup PkgPointSet2RangeSearch

computes handles to all vertices contained in the closure of the iso-rectangle `(a,b,c,d)`.

\pre `a` is the upper left point, `b` the lower left, `c` the lower right and `d` the upper right point of the iso rectangle.
The computed vertex handles will be placed as a sequence of objects in a container of value type
of `res`
which points to the first object in the sequence. The function
returns an output iterator pointing to the position beyond the end
of the sequence. `delau` is the \cgal Delaunay triangulation on which we perform the range search operation.

\cgalHeading{Requirements}

`Dt` is a \cgal Delaunay triangulation and contains the following subset of types from the concept `PointSetTraits` and from 
the Delaunay triangulation data type: 
<UL> 
<LI>`Dt::Geom_traits` 
<LI>`Dt::Vertex_handle` 
<LI>`Dt::Vertex` 
<LI>`Dt::Vertex_circulator` 
<LI>`Dt::Vertex_iterator` 
<LI>`Dt::Point` 
<LI>`Dt::Geom_traits::Circle_2` 
<LI>`Dt::Geom_traits::Bounded_side_2` 
<LI>`Dt::Geom_traits::Construct_center_2` 
<LI>`Dt::Geom_traits::Orientation_2` 
<LI>`Dt::Geom_traits::Construct_circle_2` 
</UL> 

*/
template<class Dt, class OutputIterator>
OutputIterator range_search(Dt& delau, const Dt::Point& a, const Dt::Point& b, const Dt::Point& c,const Dt::Point& d,OutputIterator res) ;

/*!
\ingroup PkgPointSet2RangeSearch

computes handles to all vertices contained in the closure of disk `C`.
The computed vertex handles will be placed as a sequence of objects in a container of value type
of `res`
which points to the first object in the sequence. The function
returns an output iterator pointing to the position beyond the end
of the sequence.
`delau` is the \cgal Delaunay triangulation on that we perform the range search operation.
` pred` controls the search operation. If `return_if_succeded` is `true`, we will end the search
after the first success of the predicate, otherwise we will continue till the search is finished.

\cgalHeading{Requirements}

For the requirements of `Dt` see the description for the non-predicate version. 

Requirements of `Pred`: 
<UL> 
<LI>`void set_result(bool);` 
<LI>`bool operator()(const Point&);` 
</UL> 
The `operator()` is used for testing the current point in the search operation. 
If this operator returns `true` and `return_if_succeded` is `true`, the range search will stop. 
Otherwise the range search operation will continue. Member function `set_result()` can be used to 
store the result of the range search in the function object. The result will be `true` if the last 
call to the `operator()` of the predicate returned `true`, `false` otherwise. 

*/
template<class Dt, class Circle, class OutputIterator, class Pred>
OutputIterator range_search(Dt& delau, const Circle& C, OutputIterator res,
Pred& pred, bool return_if_succeded);

/*!
\ingroup PkgPointSet2RangeSearch

computes handles to all vertices contained in the closure of the triangle `(a,b,c)`.

\pre `a`, `b`, and `c` must not be collinear.
The computed vertex handles will be placed as a sequence of objects in a container of value type
of `res`
which points to the first object in the sequence. The function
returns an output iterator pointing to the position beyond the end
of the sequence. 
`delau` is the \cgal Delaunay triangulation on which we perform the range search operation.
`pred` controls the search operation. If `return_if_succeded` is `true`, we will end the search
after the first success of the predicate, otherwise we will continue till the search is finished.

\cgalHeading{Requirements}

For the requirements of `Dt` see the description for the non-predicate version. 

For the requirements of `Pred` see the description above. 

*/
template<class Dt, class OutputIterator, class Pred>
OutputIterator range_search(Dt& delau, const Dt::Point& a, const Dt::Point& b, 
const Dt::Point& c,OutputIterator res, Pred& pred, bool return_if_succeded);

/*!
\ingroup PkgPointSet2RangeSearch

computes handles to all vertices contained in the closure of the iso-rectangle `(a,b,c,d)`.

\pre `a` is the upper left point, `b` the lower left, `c` the lower right and `d` the upper right point of the iso rectangle.
The computed vertex handles will be placed as a sequence of objects in a container of value type
of `res`
which points to the first object in the sequence. The function
returns an output iterator pointing to the position beyond the end
of the sequence. `delau` is the \cgal Delaunay triangulation on which we perform the range search operation.
`pred` controls the search operation. If `return_if_succeded` is `true`, we will end the search
after the first success of the predicate, otherwise we will continue till the search is finished.

\cgalHeading{Requirements}

For the requirements of `Dt` see the description for the non-predicate version. 

For the requirements of `Pred` see the description above. 

*/
template<class Dt, class OutputIterator, class Pred>
OutputIterator range_search(Dt& delau, const Dt::Point& a, const Dt::Point& b, const Dt::Point& c,const Dt::Point& d,
OutputIterator res, Pred& pred, bool return_if_succeded);

} /* namespace CGAL */
