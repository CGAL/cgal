namespace CGAL {

/*!
  \addtogroup  PkgBoxIntersectionD_box_intersection_all_pairs_d

  The function `box_intersection_all_pairs_d()` computes the pairwise intersecting boxes 
  between two sequences of iso-oriented boxes in arbitrary dimension. 
  It does so by comparing all possible pairs of boxes and is thus 
  inferior to the fast `box_intersection_d()` algorithm. 

  The sequences of boxes are given with two forward iterator ranges. The 
  sequences are not modified. For each intersecting pair of boxes a 
  `callback` function object is called with the two intersecting 
  boxes as argument; the first argument is a box from the first 
  sequence, the second argument a box from the second sequence. 

  The algorithm is interface compatible with the 
  `box_intersection_d()` function. Similarly, we call the 
  `value_type` of the iterators the <I>box handle</I>, which is 
  either our box type or a pointer type to our box type. 

  A \f$ d\f$-dimensional iso-oriented box is defined as the 
  %Cartesian product of \f$ d\f$ intervals. We call the 
  box <I>half-open</I> if the \f$ d\f$ intervals \f$ \{ [lo_i,hi_i) \,|\, 0 \leq 
  i < d\}\f$ are half-open intervals, and we call the box <I>closed</I> if 
  the \f$ d\f$ intervals \f$ \{ [lo_i,hi_i] \,|\, 0 \leq i < d\}\f$ are closed 
  intervals. Note that closed boxes support zero-width boxes and they 
  can intersect at their boundaries, while non-empty half-open boxes 
  always have a positive volume and they only intersect iff their 
  interiors overlap. The distinction between closed or half-open boxes 
  does not require a different representation of boxes, just a different 
  interpretation when comparing boxes, which is selected with the 
  `topology` parameter and its two values, 
  `Box_intersection_d::HALF_OPEN` and 
  `Box_intersection_d::CLOSED`. 

  In addition, a box has a unique `id`-number. Boxes with equal 
  `id`-number are not reported since they obviously intersect trivially. 

  The algorithm uses a traits class of the `BoxIntersectionTraits_d` 
  concept to access the boxes. A default traits class is provided that 
  assumes that the box type is a model of the `BoxIntersectionBox_d` 
  concept and that the box handle, i.e., the iterators value type, is 
  identical to the box type or a pointer to the box type. 

  An important special application of this algorithm is the test for 
  self-intersections where the second box sequence is an identical copy 
  of the first sequence including the preserved `id`-number. We 
  offer a specialized implementation 
  `box_self_intersection_all_pairs` for this application. 

  \sa \link PkgBoxIntersectionD_box_intersection_d `CGAL::box_intersection_d()` \endlink
  \sa \link PkgBoxIntersectionD_box_self_intersection_d `CGAL::box_self_intersection_d()` \endlink
  \sa \link PkgBoxIntersectionD_box_self_intersection_all_pairs_d `CGAL::box_self_intersection_all_pairs_d()` \endlink
  \sa `CGAL::Box_intersection_d::Box_traits_d<BoxHandle>` 
  \sa `BoxIntersectionBox_d` 
  \sa `BoxIntersectionTraits_d` 

\cgalHeading{Implementation}

  The algorithm is trivially testing all pairs and runs therefore in time 
  \f$ O(nm)\f$ where \f$ n\f$ is the size of the first sequence and \f$ m\f$ is the 
  size of the second sequence. 
*/

/*!
  \addtogroup  PkgBoxIntersectionD_box_intersection_d

  The function `box_intersection_d()` computes the pairwise intersecting boxes 
  between two sequences of iso-oriented boxes in arbitrary dimension. 
  The sequences of boxes are given with two random-access iterator 
  ranges and will be reordered in the course of the algorithm. For each 
  intersecting pair of boxes a `callback` function object is called 
  with the two intersecting boxes as argument; the first argument is a 
  box from the first sequence, the second argument a box from the second 
  sequence. The performance of the algorithm can be tuned with a 
  `cutoff` parameter, see the implementation section below for more 
  details. 

  The algorithm reorders the boxes in the course of the algorithm. Now, 
  depending on the size of a box it can be faster to copy the boxes, or 
  to work with pointers to boxes and copy only pointers. We offer 
  automatic support for both options. To simplify the description, let 
  us call the `value_type` of the iterators <I>box handle</I>. The 
  <I>box handle</I> can either be our box type itself or a pointer (or 
  const pointer) to the box type. 

  A \f$ d\f$-dimensional iso-oriented box is defined as the 
  %Cartesian product of \f$ d\f$ intervals. We call the 
  box <I>half-open</I> if the \f$ d\f$ intervals \f$ \{ [lo_i,hi_i) \,|\, 0 \leq 
  i < d\}\f$ are half-open intervals, and we call the box <I>closed</I> if 
  the \f$ d\f$ intervals \f$ \{ [lo_i,hi_i] \,|\, 0 \leq i < d\}\f$ are closed 
  intervals. Note that closed boxes support zero-width boxes and they 
  can intersect at their boundaries, while non-empty half-open boxes 
  always have a positive volume and they only intersect iff their 
  interiors overlap. The distinction between closed or half-open boxes 
  does not require a different representation of boxes, just a different 
  interpretation when comparing boxes, which is selected with the 
  `topology` parameter and its two values, 
  `Box_intersection_d::HALF_OPEN` and 
  `Box_intersection_d::CLOSED`. 

  In addition, a box has a unique `id`-number. It is used to order 
  boxes consistently in each dimension even if boxes have identical 
  coordinates. In consequence, the algorithm guarantees that a pair of 
  intersecting boxes is reported only once. Boxes with equal 
  `id`-number are not reported since they obviously intersect trivially. 

  The algorithm uses a traits class of the `BoxIntersectionTraits_d` 
  concept to access the boxes. A default traits class is provided that 
  assumes that the box type is a model of the `BoxIntersectionBox_d` 
  concept and that the box handle, i.e., the iterators value type, is 
  identical to the box type or a pointer to the box type. 

  An important special application of this algorithm is the test for 
  self-intersections where the second box sequence is an identical copy 
  of the first sequence including the preserved `id`-number. Note 
  that this implies that the address of the box is not sufficient for 
  the `id`-number if boxes are copied by value. To ease the use of 
  this special case we offer a simplified version of this function with 
  one iterator range only, which then creates internally the second copy 
  of the boxes, under the name `box_self_intersection_d()`. 

  In the general case, we distinguish between the bipartite case (the 
  boxes are from different sequences) and the complete case (the boxes 
  are from the same sequence, i.e., the self intersection case). The 
  default is the bipartite case, since the complete case is typically 
  handled with the simplified function call mentioned above. However, 
  the general function call offers the `setting` parameter with the 
  values `Box_intersection_d::COMPLETE` and 
  `Box_intersection_d::BIPARTITE`. 

\cgalHeading{Requirements}

  <UL> 
  <LI>`RandomAccessIterator1`, and \f$ \ldots\f$ `2`, must be 
  mutable random-access iterators and both value types must be 
  the same. We call this value type `Box_handle` in the following. 
  <LI>`Callback` must be of the `BinaryFunction` concept. 
  The `Box_handle` must be convertible to both argument types. The 
  return type is not used and can be `void`. 
  <LI>The `Box_handle` must be a model of the `Assignable` concept. 
  <LI>In addition, if the default box traits is used the `Box_handle` must 
  be a class type `T` or a pointer to a class type `T`, where 
  `T` must be a model of the `BoxIntersectionBox_d` concept. 
  In both cases, the default box traits specializes to a suitable 
  implementation. 
  <LI>`BoxTraits` must be of the `BoxIntersectionTraits_d` concept. 
  </UL> 

  \sa \link PkgBoxIntersectionD_box_self_intersection_d `CGAL::box_self_intersection_d()` \endlink
  \sa \link PkgBoxIntersectionD_box_intersection_all_pairs_d `CGAL::box_intersection_all_pairs_d()` \endlink
  \sa `CGAL::Box_intersection_d::Box_traits_d<BoxHandle>` 
  \sa `BoxIntersectionBox_d` 
  \sa `BoxIntersectionTraits_d` 

\cgalHeading{Implementation}

  The implemented algorithm is described in \cgalCite{cgal:ze-fsbi-02} as 
  version two. Its performance depends on a `cutoff` parameter. 
  When the size of both iterator ranges drops below the `cutoff` 
  parameter the function switches from the streamed segment-tree 
  algorithm to the two-way-scan algorithm, see \cgalCite{cgal:ze-fsbi-02} 
  for the details. 

  The streamed segment-tree algorithm needs \f$ O(n \log^d (n) + k)\f$ 
  worst-case running time and \f$ O(n)\f$ space, where \f$ n\f$ is the number of 
  boxes in both input sequences, \f$ d\f$ the (constant) dimension of the 
  boxes, and \f$ k\f$ the output complexity, i.e., the number of pairwise 
  intersections of the boxes. The two-way-scan algorithm needs \f$ O(n \log 
  (n) + l)\f$ worst-case running time and \f$ O(n)\f$ space, where \f$ l\f$ is the 
  number of pairwise overlapping intervals in one dimensions (the 
  dimension where the algorithm is used instead of the segment tree). 
  Note that \f$ l\f$ is not necessarily related to \f$ k\f$ and using the 
  two-way-scan algorithm is a heuristic. 

  Unfortunately, we have no general method to automatically determine an 
  optimal cutoff parameter, since it depends on the used hardware, the 
  runtime ratio between callback runtime and segment-tree runtime, and 
  of course the number of boxes to be checked and their distribution. In 
  cases where the callback runtime is dominant, it may be best to make 
  the threshold parameter small. Otherwise a `cutoff`\f$ =\sqrt{n}\f$ can 
  lead to acceptable results. For well distributed boxes the original 
  paper \cgalCite{cgal:ze-fsbi-02} gives optimal cutoffs in the thousands. 
  Anyway, for optimal runtime some experiments to compare different 
  cutoff parameters are recommended. See also 
  Section \ref secboxintersperformance . 

\cgalHeading{Example}

  The box implementation provided with 
  `Box_intersection_d::Box_d<double,2>` has a special 
  constructor for the \cgal bounding box type `Bbox_2` (and 
  similar for dimension 3). We use this in the example to create \f$ 3 
  \times 3\f$ `boxes` in a grid layout. Additionally we pick the 
  center box and the box in the upper-right corner as our second box 
  sequence `query`. 

  The default policy of the box type implements the `id`-number with 
  an explicit counter in the boxes, which is the default choice since it 
  always works. We use the `id`-number in our callback function to 
  report the result of the intersection algorithm call. The result will 
  be that the first `query` box intersects all nine `boxes` and 
  the second `query` box intersects the four boxes in the 
  upper-right quadrant. 

  \cgalExample{Box_intersection_d/minimal.cpp} 

*/

/*!
  \ingroup PkgBoxIntersectionD_box_intersection_all_pairs_d

  Invocation of box intersection with default box traits
  `Box_intersection_d::Box_traits_d<Box_handle>`, where
  `Box_handle` corresponds to the iterator value type of
  `ForwardIterator1`.

*/
template< class ForwardIterator1, 
          class ForwardIterator2, 
          class Callback >
void box_intersection_all_pairs_d(
  ForwardIterator1 begin1, ForwardIterator1 end1,
  ForwardIterator2 begin2, ForwardIterator2 end2,
  Callback callback,
  CGAL::Box_intersection_d::Topology topology = CGAL::Box_intersection_d::CLOSED);

/*!
  \ingroup PkgBoxIntersectionD_box_intersection_all_pairs_d

  Invocation with custom box traits.

*/
template< class ForwardIterator1,
          class ForwardIterator2,
          class Callback, class BoxTraits >
void box_intersection_all_pairs_d(
  ForwardIterator1 begin1, ForwardIterator1 end1,
  ForwardIterator2 begin2, ForwardIterator2 end2,
  Callback callback,
  BoxTraits box_traits,
  CGAL::Box_intersection_d::Topology topology = CGAL::Box_intersection_d::CLOSED);

} /* namespace CGAL */

namespace CGAL {

/*!
  \ingroup PkgBoxIntersectionD_box_intersection_d

  Invocation of box intersection with default box traits
  `Box_intersection_d::Box_traits_d<Box_handle>`, where
  `Box_handle` corresponds to the iterator value type of
  `RandomAccessIterator1`.


*/
template< class RandomAccessIterator1, 
          class RandomAccessIterator2, 
          class Callback >
void box_intersection_d(
  RandomAccessIterator1 begin1, RandomAccessIterator1 end1,
  RandomAccessIterator2 begin2, RandomAccessIterator2 end2,
  Callback callback,
  std::ptrdiff_t cutoff = 10,
  CGAL::Box_intersection_d::Topology topology = CGAL::Box_intersection_d::CLOSED,
  CGAL::Box_intersection_d::Setting setting = CGAL::Box_intersection_d::BIPARTITE);

/*!
  \ingroup PkgBoxIntersectionD_box_intersection_d

  Invocation with custom box traits.

*/
template< class RandomAccessIterator1,
          class RandomAccessIterator2,
          class Callback, class BoxTraits >
void box_intersection_d(
  RandomAccessIterator1 begin1, RandomAccessIterator1 end1,
  RandomAccessIterator2 begin2, RandomAccessIterator2 end2,
  Callback callback,
  BoxTraits box_traits,
  std::ptrdiff_t cutoff = 10,
  CGAL::Box_intersection_d::Topology topology = CGAL::Box_intersection_d::CLOSED,
  CGAL::Box_intersection_d::Setting setting = CGAL::Box_intersection_d::BIPARTITE);

} /* namespace CGAL */

namespace CGAL {

/*!
  \addtogroup  PkgBoxIntersectionD_box_self_intersection_all_pairs_d
  The function `box_self_intersection_all_pairs_d()` computes the pairwise intersecting boxes 
  in a sequence of iso-oriented boxes in arbitrary dimension. 
  It does so by comparing all possible pairs of boxes and is thus 
  inferior to the fast `box_self_intersection_d()` algorithm. 

  The sequence of boxes is given with a forward iterator range. The 
  sequences are not modified. For each intersecting pair of boxes a 
  `callback` function object is called with the two intersecting 
  boxes as argument. 

  The algorithm is interface compatible with the 
  `box_self_intersection_d()` function. Similarly, we call the 
  `value_type` of the iterators the <I>box handle</I>, which is 
  either our box type or a pointer type to our box type. 

  A \f$ d\f$-dimensional iso-oriented box is defined as the 
  %Cartesian product of \f$ d\f$ intervals. We call the 
  box <I>half-open</I> if the \f$ d\f$ intervals \f$ \{ [lo_i,hi_i) \,|\, 0 \leq 
  i < d\}\f$ are half-open intervals, and we call the box <I>closed</I> if 
  the \f$ d\f$ intervals \f$ \{ [lo_i,hi_i] \,|\, 0 \leq i < d\}\f$ are closed 
  intervals. Note that closed boxes support zero-width boxes and they 
  can intersect at their boundaries, while non-empty half-open boxes 
  always have a positive volume and they only intersect iff their 
  interiors overlap. The distinction between closed or half-open boxes 
  does not require a different representation of boxes, just a different 
  interpretation when comparing boxes, which is selected with the 
  `topology` parameter and its two values, 
  `Box_intersection_d::HALF_OPEN` and 
  `Box_intersection_d::CLOSED`. 

  The algorithm uses a traits class of the `BoxIntersectionTraits_d` 
  concept to access the boxes. A default traits class is provided that 
  assumes that the box type is a model of the `BoxIntersectionBox_d` 
  concept and that the box handle, i.e., the iterators value type, is 
  identical to the box type or a pointer to the box type. 

\cgalHeading{Requirements}

  <UL> 
  <LI>`ForwardIterator` must be a forward iterator. We call its 
  value type `Box_handle` in the following. 
  <LI>`Callback` must be of the `BinaryFunction` concept. 
  The `Box_handle` must be convertible to both argument types. The 
  return type is not used and can be `void`. 
  <LI>The `Box_handle` must be a model of the `Assignable` concept. 
  <LI>In addition, if the default box traits is used the `Box_handle` must 
  be a class type `T` or a pointer to a class type `T`, where 
  `T` must be a model of the `BoxIntersectionBox_d` concept. 
  In both cases, the default box traits specializes to a suitable 
  implementation. 
  <LI>`BoxTraits` must be of the `BoxIntersectionTraits_d` concept. 
  </UL> 

  \sa \link PkgBoxIntersectionD_box_intersection_d `CGAL::box_intersection_d()` \endlink
  \sa \link PkgBoxIntersectionD_box_self_intersection_d `CGAL::box_self_intersection_d()` \endlink
  \sa \link PkgBoxIntersectionD_box_intersection_all_pairs_d `CGAL::box_intersection_all_pairs_d()` \endlink
  \sa `CGAL::Box_intersection_d::Box_traits_d<BoxHandle>` 
  \sa `BoxIntersectionBox_d` 
  \sa `BoxIntersectionTraits_d` 

\cgalHeading{Implementation}

  The algorithm is trivially testing all pairs and runs therefore in time 
  \f$ O(n^2)\f$ where \f$ n\f$ is the size of the input sequence. This algorithm 
  does not use the id-number of the boxes. 

*/

/*!
  \ingroup PkgBoxIntersectionD_box_self_intersection_all_pairs_d

  Invocation of box intersection with default box traits
  `Box_intersection_d::Box_traits_d<Box_handle>`, where
  `Box_handle` corresponds to the iterator value type of
  `ForwardIterator`.


*/
template< class ForwardIterator, class Callback >
void box_self_intersection_all_pairs_d(
  ForwardIterator begin, ForwardIterator end,
  Callback callback,
  CGAL::Box_intersection_d::Topology topology = CGAL::Box_intersection_d::CLOSED);

/*!
  \ingroup PkgBoxIntersectionD_box_self_intersection_all_pairs_d
  Invocation with custom box traits.


*/
template< class ForwardIterator,
          class Callback, class BoxTraits >
void box_self_intersection_all_pairs_d(
  ForwardIterator begin, ForwardIterator end,
  Callback callback,
  BoxTraits box_traits,
  CGAL::Box_intersection_d::Topology topology = CGAL::Box_intersection_d::CLOSED);

} /* namespace CGAL */

namespace CGAL {

/*!
  \addtogroup  PkgBoxIntersectionD_box_self_intersection_d 
  The function `box_self_intersection_d()` computes the pairwise intersecting boxes 
  in a sequence of iso-oriented boxes in arbitrary dimension. 
  The sequence of boxes is given with as a random-access iterator 
  range and will be reordered in the course of the algorithm. For each 
  intersecting pair of boxes a `callback` function object is called 
  with the two intersecting boxes as argument; the first argument is a 
  box from the sequence, the second argument is a copy of a box from the 
  sequence. The performance of the algorithm can be tuned with a 
  `cutoff` parameter, see the implementation section of the 
  `box_intersection_d()` function. 

  The algorithm creates a second copy of the boxes and reorders the 
  boxes in the course of the algorithm. Now, depending on the size of a 
  box it can be faster to copy the boxes, or to work with pointers to 
  boxes and copy only pointers. We offer automatic support for both 
  options. To simplify the description, let us call the `value_type` 
  of the iterators <I>box handle</I>. The <I>box handle</I> can 
  either be our box type itself or a pointer (or const pointer) to the 
  box type. 

  A \f$ d\f$-dimensional iso-oriented box is defined as the 
  %Cartesian product of \f$ d\f$ intervals. We call the 
  box <I>half-open</I> if the \f$ d\f$ intervals \f$ \{ [lo_i,hi_i) \,|\, 0 \leq 
  i < d\}\f$ are half-open intervals, and we call the box <I>closed</I> if 
  the \f$ d\f$ intervals \f$ \{ [lo_i,hi_i] \,|\, 0 \leq i < d\}\f$ are closed 
  intervals. Note that closed boxes support zero-width boxes and they 
  can intersect at their boundaries, while non-empty half-open boxes 
  always have a positive volume and they only intersect iff their 
  interiors overlap. The distinction between closed or half-open boxes 
  does not require a different representation of boxes, just a different 
  interpretation when comparing boxes, which is selected with the 
  `topology` parameter and its two values, 
  `Box_intersection_d::HALF_OPEN` and 
  `Box_intersection_d::CLOSED`. 

  In addition, a box has a unique `id`-number. It is used to order 
  boxes consistently in each dimension even if boxes have identical 
  coordinates. In consequence, the algorithm guarantees that a pair of 
  intersecting boxes is reported only once. This self-intersection 
  function creates internally a second copy of the box sequence. The 
  copying has to preserve the `id`-number of boxes. Note that this 
  implies that the address of the box is not sufficient for the 
  `id`-number if boxes are copied by value. Boxes of equal 
  `id`-number are not reported as intersecting pairs since they are 
  always intersecting trivially. 

  The algorithm uses a traits class of the `BoxIntersectionTraits_d` 
  concept to access the boxes. A default traits class is provided that 
  assumes that the box type is a model of the `BoxIntersectionBox_d` 
  concept and that the box handle, i.e., the iterators value type, is 
  identical to the box type or a pointer to the box type. 

\cgalHeading{Requirements}

  <UL> 
  <LI>`RandomAccessIterator` must be a mutable random-access 
  iterator. We call its value type `Box_handle` in the following. 
  <LI>`Callback` must be of the `BinaryFunction` concept. 
  The `Box_handle` must be convertible to both argument types. The 
  return type is not used and can be `void`. 
  <LI>The `Box_handle` must be a model of the `Assignable` concept. 
  <LI>In addition, if the default box traits is used the `Box_handle` must 
  be a class type `T` or a pointer to a class type `T`, where 
  `T` must be a model of the `BoxIntersectionBox_d` concept. 
  In both cases, the default box traits specializes to a suitable 
  implementation. 
  <LI>`BoxTraits` must be of the `BoxIntersectionTraits_d` concept. 
  </UL> 

  \sa \link PkgBoxIntersectionD_box_intersection_d `CGAL::box_intersection_d()` \endlink
  \sa \link PkgBoxIntersectionD_box_self_intersection_all_pairs_d `CGAL::box_self_intersection_all_pairs_d()` \endlink
  \sa `CGAL::Box_intersection_d::Box_traits_d<BoxHandle>` 
  \sa `BoxIntersectionBox_d` 
  \sa `BoxIntersectionTraits_d` 

\cgalHeading{Implementation}

  See the implementation section of the `box_intersection_d()` 
  function.

\cgalHeading{Example}

  The box implementation provided with 
  `Box_intersection_d::Box_d<double,2>` has a special 
  constructor for the \cgal bounding box type `Bbox_2` (and 
  similar for dimension 3). We use this in the example to create \f$ 3 
  \times 3\f$ `boxes` in a grid layout. 

  The default policy of the box type implements the `id`-number with 
  an explicit counter in the boxes, which is the default choice since it 
  always works. We use the `id`-number in our callback function to 
  report the result of the intersection algorithm call. The result will 
  be 20 pairwise intersections, but the order in which they are reported 
  is non-intuitive. 

  \cgalExample{Box_intersection_d/minimal_self.cpp} 

*/

/*!
  \ingroup PkgBoxIntersectionD_box_self_intersection_d

  Invocation of box intersection with default box traits
  `Box_intersection_d::Box_traits_d<Box_handle>`, where
  `Box_handle` corresponds to the iterator value type of
  `RandomAccessIterator`.

*/
template< class RandomAccessIterator, class Callback >
void box_self_intersection_d(
  RandomAccessIterator begin, RandomAccessIterator end,
  Callback callback,
  std::ptrdiff_t cutoff = 10,
  CGAL::Box_intersection_d::Topology topology = CGAL::Box_intersection_d::CLOSED);

/*!
  \ingroup PkgBoxIntersectionD_box_self_intersection_d

  Invocation with custom box traits.
*/
template< class RandomAccessIterator,
          class Callback, class BoxTraits >
void box_self_intersection_d(
  RandomAccessIterator begin, RandomAccessIterator end,
  Callback callback,
  BoxTraits box_traits,
  std::ptrdiff_t cutoff = 10,
  CGAL::Box_intersection_d::Topology topology = CGAL::Box_intersection_d::CLOSED);

} /* namespace CGAL */
