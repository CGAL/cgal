namespace CGAL {

/*!
\addtogroup boolean_complement Complement Functions
\ingroup PkgBooleanSetOperations2
\anchor ref_bso_complement 

The `complement` function is overloaded. Depending on the
type of polygon `pgn` the complement is either a single (general) polygon with 
holes, or several (general) poylgons with holes. In the latter case 
the `complement function` writes them into an output iterator 
`oi`.

\param pgn The input polygon for the `complement` function. It may be of the type
`Polygon_2`, `General_polygon_2`, `Polygon_with_holes_2`, or
`General_polygon_with_holes_2`.



\sa \link boolean_do_intersect `CGAL::do_intersect()` \endlink
\sa \link boolean_intersection `CGAL::intersection()` \endlink
\sa \link boolean_join `CGAL::join()` \endlink
\sa \link boolean_difference `CGAL::difference()` \endlink
\sa \link boolean_symmetric_difference `CGAL::symmetric_difference()` \endlink
*/
/// @{

/*!
 writes the complement of the polygon `pgn` into the  polygon with holes `res`.
 */
template <class Kernel, class Container>
void complement(const Polygon_2<Kernel, Container> & pgn, Polygon_with_holes_2<Kernel, Container> & res);

/*!
  writes the complement of the general polygon `pgn` into the  general polygon with holes `res`.
 */
template <class Traits>
void complement(const General_polygon_2<Traits> & pgn, General_polygon_with_holes_2<Traits> & res);

/*!
  writes the complement of the polygon with holes `pgn` into the  output iterator `oi`.
  The value type of `oi` is `Polygon_with_holes_2`.
 */
template <class Traits, class OutputIterator>
OutputIterator complement(const Polygon_with_holes_2<Kernel, Container> & pgn, OutputIterator oi);

/*!
  writes the complement of the general polygon with holes `pgn` into the  output iterator `oi`.
  The value type of `oi` is `General_polygon_with_holes_2`.
 */
template <class Traits, class OutputIterator>
OutputIterator complement(const General_polygon_with_holes_2<General_polygon_2<Traits> > & pgn, OutputIterator oi);
/// @}

} /* namespace CGAL */

namespace CGAL {

/*!
\addtogroup boolean_difference Difference Functions
\ingroup PkgBooleanSetOperations2
\anchor ref_bso_difference 

Each one of these functions computes the difference between two given 
polygons `p1` and `p2`, and inserts the resulting polygons 
with holes into an output container through the output iterator `oi`. 
The value type of the `OutputIterator` is either 
`Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

The signature of the function is
\code
  OutputIterator difference(const Type1 & p1, const Type2 & p2, OutputIterator oi);
\endcode

\cgalHeading{Parameters}

The types of the paramters of the `difference()` function are any of the following combinations.

<div align="left"> 
<table cellpadding=3 border="1"> 
<tr><th> Type1</th><th>Type2</th></tr> 
<tr><td valign="center">Polygon_2</td><td valign="center">Polygon_2</td></tr> 
<tr><td valign="center">Polygon_2</td><td valign="center">Polygon_with_holes_2</td></tr> 
<tr><td valign="center">Polygon_with_holes_2</td><td valign="center">Polygon_2</td></tr> 
<tr><td valign="center">Polygon_with_holes_2</td><td valign="center">Polygon_with_holes_2</td></tr> 
<tr><td valign="center">General_polygon_2</td><td valign="center">General_polygon_2</td></tr> 
<tr><td valign="center">General_polygon_2</td><td valign="center">General_polygon_with_holes_2</td></tr> 
<tr><td valign="center">General_polygon_with_holes_2</td><td valign="center">General_polygon_2</td></tr> 
<tr><td valign="center">General_polygon_with_holes_2</td><td valign="center">General_polygon_with_holes_2</td></tr> 
</table> 
</div> 

\sa \link boolean_do_intersect `CGAL::do_intersect()` \endlink
\sa \link boolean_intersection `CGAL::intersection()` \endlink
\sa \link boolean_join `CGAL::join()` \endlink
\sa \link boolean_symmetric_difference `CGAL::symmetric_difference()` \endlink

*/


/// @{
/*!
 writes the difference of the polygons `p1` and `p2` into the output iterator `oi`.
 The value type of `oi` is `Polygon_with_holes_2`.
 */
template <class Kernel, class Container, class OutputIterator>
OutputIterator difference(const Polygon_2<Kernel, Container> & p1,
                          const Polygon_2<Kernel, Container> & p2,
                          OutputIterator oi);

/*!
 writes the difference of the polygons `p1` and `p2` into the output iterator `oi`.
 The value type of `oi` is `Polygon_with_holes_2`.
 */
template <class Kernel, class Container, class OutputIterator>
OutputIterator difference(const Polygon_2<Kernel, Container> & p1,
                          const Polygon_with_holes_2<Kernel,Container> & p2,
                          OutputIterator oi);

/*!
 writes the difference of the polygons `p1` and `p2` into the output iterator `oi`.
 The value type of `oi` is `Polygon_with_holes_2`.
 */
template <class Kernel, class Container, class OutputIterator>
OutputIterator difference(const Polygon_with_holes_2<Kernel, Container> & p1,
                          const Polygon_2<Kernel, Container> & p2,
                          OutputIterator oi);

/*!
 writes the difference of the polygons `p1` and `p2` into the output iterator `oi`.
 The value type of `oi` is `Polygon_with_holes_2`.
 */
template <class Kernel, class Container, class OutputIterator>
OutputIterator difference(const Polygon_with_holes_2<Kernel, Container> & p1,
                          const Polygon_with_holes_2<Kernel, Container> & p2,
                          OutputIterator oi);

/*!
 writes the difference of the general polygons `p1` and `p2` into the output iterator `oi`.
 The value type of `oi` is `General_polygon_with_holes_2`.
 */
template <class Traits, class OutputIterator>
OutputIterator difference(const General_polygon_2<Traits> & p1,
                          const General_polygon_2<Traits> & p2,
                          OutputIterator oi);

/*!
 writes the difference of the general polygons `p1` and `p2` into the output iterator `oi`.
 The value type of `oi` is `General_polygon_with_holes_2`.
 */
template <class Traits, class OutputIterator>
OutputIterator difference(const General_polygon_with_holes_2<General_polygon_2<Traits> > & p1,
                          const General_polygon_2<Traits> & p2,
                          OutputIterator oi);


/*!
 writes the difference of the general polygons `p1` and `p2` into the output iterator `oi`.
 The value type of `oi` is `General_polygon_with_holes_2`.
 */
template <class Traits, class OutputIterator>
OutputIterator difference(const General_polygon_2<Traits> & p1,
                          const General_polygon_with_holes_2<General_polygon_2<Traits> > & p2,
                          OutputIterator oi);


/*!
 writes the difference of the general polygons `p1` and `p2` into the output iterator `oi`.
 The value type of `oi` is `General_polygon_with_holes_2`.
 */
template <class Polygon, class OutputIterator>
OutputIterator difference(const General_polygon_with_holes_2<Polygon> & p1,
                          const General_polygon_with_holes_2<Polygon> & p2,
                          OutputIterator oi);
/// @}

} /* namespace CGAL */

namespace CGAL {

/*!
\addtogroup boolean_do_intersect Intersection Testing Functions
\ingroup PkgBooleanSetOperations2
\anchor ref_bso_do_intersect 

Each one of these functions computes if the interior of two given 
polygons `p1` and `p2` intersect. 

The signature of the function is
\code
  bool do_intersect(const Type1 & p1, const Type2 & p2);
\endcode

\cgalHeading{Parameters}

The types of the paramters of the `do_intersect()` function are any of the following combinations.

<div align="left"> 
<table cellpadding=3 border="1"> 
<tr><th> Type1 </th><th>Type2</th></tr> 
<tr><td valign="center">Polygon_2</td><td valign="center">Polygon_2</td></tr> 
<tr><td valign="center">Polygon_2</td><td valign="center">Polygon_with_holes_2</td></tr> 
<tr><td valign="center">Polygon_with_holes_2</td><td valign="center">Polygon_2</td></tr> 
<tr><td valign="center">Polygon_with_holes_2</td><td valign="center">Polygon_with_holes_2</td></tr> 
<tr><td valign="center">General_polygon_2</td><td valign="center">General_polygon_2</td></tr> 
<tr><td valign="center">General_polygon_2</td><td valign="center">General_polygon_with_holes_2</td></tr> 
<tr><td valign="center">General_polygon_with_holes_2</td><td valign="center">General_polygon_2</td></tr> 
<tr><td valign="center">General_polygon_with_holes_2</td><td valign="center">General_polygon_with_holes_2</td></tr> 
</table> 
</div> 


\sa \link boolean_intersection `CGAL::intersection()` \endlink
\sa \link boolean_join `CGAL::join()` \endlink
\sa \link boolean_difference `CGAL::difference()` \endlink
\sa \link boolean_symmetric_difference `CGAL::symmetric_difference()` \endlink

*/

/// @{

/*!
  returns `true` if the polygons `p1` and `p2` intersect in their interior.
 */
template <class Kernel, class Container>
bool do_intersect(const Polygon_2<Kernel, Container> & p1,
                  const Polygon_2<Kernel, Container> & p2);

/*!
  returns `true` if the polygons `p1` and `p2` intersect in their interior.
 */
template <class Kernel, class Container>
bool do_intersect(const Polygon_2<Kernel, Container> & p1,
                  const Polygon_with_holes_2<Kernel, Container> & p2);

/*!
  returns `true` if the polygons `p1` and `p2` intersect in their interior.
  returns `true` if the interior of polygons `p1` and `p2` intersect.
 */
template <class Kernel, class Container>
bool do_intersect(const Polygon_with_holes_2<Kernel, Container> & p1,
                  const Polygon_2<Kernel, Container> & p2);

/*!
  returns `true` if the polygons `p1` and `p2` intersect in their interior.

 */
template <class Kernel, class Container>
bool do_intersect(const Polygon_with_holes_2<Kernel, Container> & p1,
                  const Polygon_with_holes_2<Kernel, Container> & p2);

/*!
  returns `true` if the general polygons `p1` and `p2` intersect in their interior.
 */
template <class Traits>
bool do_intersect(const General_polygon_2<Traits> & p1,
                  const General_polygon_2<Traits> & p2);
  
/*!
  returns `true` if the general polygons `p1` and `p2` intersect in their interior.
 */
template <class Traits>
bool do_intersect(const General_polygon_2<Traits> & p1,
                  const General_polygon_with_holes_2<General_polygon_2<Traits> > & p2);

/*!
  returns `true` if the general polygons `p1` and `p2` intersect in their interior.
 */
template <class Traits>
bool do_intersect(const General_polygon_with_holes_2<General_polygon_2<Traits> > & p1,
                  const General_polygon_2<Traits> & p2);


/*!
  returns `true` if the general polygons `p1` and `p2` intersect in their interior.
 */
template <class Polygon>
bool do_intersect(const General_polygon_with_holes_2<Polygon> & p1,
                  const General_polygon_with_holes_2<Polygon> & p2);

  /*!
    returns `true`, if the set of general polygons (or general 
    polygons with holes) in the given range intersect in their interior, 
    and `false` otherwise. (The value type of the input iterator is 
    used to distinguish between the two). 
  */
template <class InputIterator>
bool do_intersect(InputIterator begin, InputIterator end);

  /*!
    returns `true`, if the set of general polygons and general 
    polygons with holes in the given two ranges respectively intersect in 
    their interior, and `false` otherwise. 
  */
template <class InputIterator1, class InputIterator2>
bool do_intersect(InputIterator1 pgn_begin1,
                  InputIterator1 pgn_end1,
                  InputIterator2 pgn_begin2,
                  InputIterator2 pgn_end2);
/// @}
} /* namespace CGAL */

namespace CGAL {

/*!
\addtogroup boolean_intersection Intersection Functions
\ingroup PkgBooleanSetOperations2
\anchor ref_bso_intersection

Each one of these functions computes the intersection of two given 
polygons `p1` and `p2`, inserts the resulting polygons with 
holes into an output container through a given output iterator 
`oi`, and returns the output iterator. The value type of the 
`OutputIterator` is either `Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 


The signature of the function is
\code
  OutputIterator intersection(const Type1 & p1, const Type2 & p2, OutputIterator oi);
\endcode


\cgalHeading{Parameters}

The types of the paramters of the `intersection()` function are any of the following combinations.


<div align="left"> 
<table cellpadding=3 border="1"> 
<tr><th> Type1</th><th> Type2</th></tr> 
<tr><td valign="center">Polygon_2</td><td valign="center">Polygon_2</td></tr> 
<tr><td valign="center">Polygon_2</td><td valign="center">Polygon_with_holes_2</td></tr> 
<tr><td valign="center">Polygon_with_holes_2</td><td valign="center">Polygon_2</td></tr> 
<tr><td valign="center">Polygon_with_holes_2</td><td valign="center">Polygon_with_holes_2</td></tr> 
<tr><td valign="center">General_polygon_2</td><td valign="center">General_polygon_2</td></tr> 
<tr><td valign="center">General_polygon_2</td><td valign="center">General_polygon_with_holes_2</td></tr> 
<tr><td valign="center">General_polygon_with_holes_2</td><td valign="center">General_polygon_2</td></tr> 
<tr><td valign="center">General_polygon_with_holes_2</td><td valign="center">General_polygon_with_holes_2</td></tr> 
</table> 
</div> 


\sa \link boolean_do_intersect `CGAL::do_intersect()` \endlink
\sa \link boolean_join `CGAL::join()` \endlink
\sa \link boolean_difference `CGAL::difference()` \endlink
\sa \link boolean_symmetric_difference `CGAL::symmetric_difference()` \endlink

*/


/// @{

/*!
 writes the intersection of the polygons `p1` and `p2` into the output iterator `oi`.
 The value type of `oi` is `Polygon_with_holes_2`.
 */
OutputIterator intersection(const Type1 & p1, const Type2 & p2,
                            OutputIterator oi);

/*!
 writes the intersection of the polygons `p1` and `p2` into the output iterator `oi`.
 The value type of `oi` is `Polygon_with_holes_2`.
 */
template <class Kernel, class Container, class OutputIterator>
OutputIterator intersection(const Polygon_2<Kernel, Container> & p1,
                            const Polygon_2<Kernel, Container> & p2,
                            OutputIterator oi);

/*!
 writes the intersection of the polygons `p1` and `p2` into the output iterator `oi`.
 The value type of `oi` is `Polygon_with_holes_2`.
 */
template <class Kernel, class Container, class OutputIterator>
OutputIterator intersection(const Polygon_2<Kernel, Container> & p1,
                            const Polygon_with_holes_2<Kernel, Container> & p2,
                            OutputIterator oi);

/*!
 writes the intersection of the polygons `p1` and `p2` into the output iterator `oi`.
 The value type of `oi` is `Polygon_with_holes_2`.
 */
template <class Kernel, class Container, class OutputIterator>
OutputIterator intersection(const Polygon_with_holes_2<Kernel, Container> & p1,
                            const Polygon_2<Kernel, Container> & p2,
                            OutputIterator oi);

/*!
 writes the intersection of the polygons `p1` and `p2` into the output iterator `oi`.
 The value type of `oi` is `Polygon_with_holes_2`.
 */
template <class Kernel, class Container, class OutputIterator>
OutputIterator intersection(const Polygon_with_holes_2<Kernel, Container> & p1,
                            const Polygon_with_holes_2<Kernel, Container> & p2,
                            OutputIterator oi);

/*!
 writes the intersection of the general polygons `p1` and `p2` into the output iterator `oi`.
 The value type of `oi` is `General_polygon_with_holes_2`.
 */
template <class Traits, class OutputIterator>
OutputIterator intersection(const General_polygon_2<Traits> & p1,
const General_polygon_2<Traits> & p2,
OutputIterator oi);

/*!
 writes the intersection of the general polygons `p1` and `p2` into the output iterator `oi`.
 The value type of `oi` is `General_polygon_with_holes_2`.
 */
template <class Traits, class OutputIterator>
OutputIterator intersection(const General_polygon_with_holes_2<General_polygon_2<Traits> > & p1,
                            const General_polygon_2<Traits> & p2,
                            OutputIterator oi);

/*!
 writes the intersection of the general polygons `p1` and `p2` into the output iterator `oi`.
 The value type of `oi` is `General_polygon_with_holes_2`.
 */
template <class Traits, class OutputIterator>
OutputIterator intersection(const General_polygon_2<Traits> & p1,
                            const General_polygon_with_holes_2<General_polygon_2<Traits> > & p2,
                            OutputIterator oi);

/*!
 writes the intersection of the general polygons `p1` and `p2` into the output iterator `oi`.
 The value type of `oi` is `General_polygon_with_holes_2`.
 */
template <class Polygon, class OutputIterator>
OutputIterator intersection(const General_polygon_with_holes_2<Polygon> & p1,
                            const General_polygon_with_holes_2<Polygon> & p2,
                            OutputIterator oi);


/*!
  computes the intersection of the general polygons (or general polygons with 
  holes) in the given range. (The value type of the input iterator is 
  used to distinguish between the two.) The result, represented by a set 
  of general polygon with holes, is written into the output iterator `oi`. 
  The output iterator is returned. The value type of the `OutputIterator` is 
  `Traits::Polygon_with_holes_2`. 
*/
template <class InputIterator, class OutputIterator>
OutputIterator intersection(InputIterator begin, InputIterator end,
                            OutputIterator oi);

/*!
  computes the intersection of the general polygons and general polygons 
  with holes in the given two ranges. The result, represented by a set 
  of general polygon with holes, is written into the output iterator `oi`. 
  The output iterator is returned. The value type of the `OutputIterator` is 
  `Traits::Polygon_with_holes_2`. 
*/
template <class InputIterator1, class InputIterator2,
class OutputIterator>
OutputIterator intersection(InputIterator1 pgn_begin1,
InputIterator1 pgn_end1,
InputIterator2 pgn_begin2,
InputIterator2 pgn_end2,
OutputIterator oi);

/// @}

} /* namespace CGAL */

namespace CGAL {

/*!
\addtogroup boolean_join Union Functions
\ingroup PkgBooleanSetOperations2
\anchor ref_bso_union

Each one of these functions computes the union of two given polygons 
`p1` and `p2`. If the two given polygons overlap, it returns 
`true`, and places the resulting polygon in `p`. Otherwise, it 
returns `false`. 


The signature of the function is
\code
  bool join(const Type1 & p1, const Type2 & p2, General_polygon_with_holes_2 & res);
\endcode

\cgalHeading{Parameters}

The types of the paramters of the `join()` function are any of the following combinations.

<div align="left"> 
<table cellpadding=3 border="1"> 
<tr><th> Type1</th><th> Type2</th></tr> 
<tr><td valign="center">Polygon_2</td><td valign="center">Polygon_2</td></tr> 
<tr><td valign="center">Polygon_2</td><td valign="center">polygon_with_holes_2</td></tr> 
<tr><td valign="center">Polygon_with_holes_2</td><td valign="center">Polygon_2</td></tr> 
<tr><td valign="center">Polygon_with_holes_2</td><td valign="center">Polygon_with_holes_2</td></tr> 
<tr><td valign="center">General_polygon_2</td><td valign="center">General_polygon_2</td></tr> 
<tr><td valign="center">General_polygon_2</td><td valign="center">General_polygon_with_holes_2</td></tr> 
<tr><td valign="center">General_polygon_with_holes_2</td><td valign="center">General_polygon_2</td></tr> 
<tr><td valign="center">General_polygon_with_holes_2</td><td valign="center">General_polygon_with_holes_2</td></tr> 
</table> 
</div> 

\sa \link boolean_do_intersect `CGAL::do_intersect()` \endlink
\sa \link boolean_intersection `CGAL::intersection()` \endlink
\sa \link boolean_difference `CGAL::difference()` \endlink
\sa \link boolean_symmetric_difference `CGAL::symmetric_difference()` \endlink

*/

/// @{

/*!
 writes the union of the polygons `p1` and `p2` into the polygon with holes `res`.
 Returns `true` if the two given polygons overlap.
 */
template <class Kernel, class Container>
bool join(const Polygon_2<Kernel, Container> & p1,
          const Polygon_2<Kernel, Container> & p2,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container> > & res);


/*!
 writes the union of the polygons `p1` and `p2` into the polygon with holes `res`.
 Returns `true` if the two given polygons overlap.
 */
template <class Kernel, class Container>
bool join(const Polygon_2<Kernel, Container> & p1,
          const Polygon_with_holes_2<Kernel,Container> & p2,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container> > & res);


/*!
 writes the union of the polygons `p1` and `p2` into the polygon with holes `res`.
 Returns `true` if the two given polygons overlap.
 */
template <class Kernel, class Container>
bool join(const Polygon_with_holes_2<Kernel, Container> & p2,
          const Polygon_2<Kernel, Container> & p1,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container> > & res);


/*!
 writes the union of the polygons `p1` and `p2` into the polygon with holes `res`.
 Returns `true` if the two given polygons overlap.
 */
template <class Kernel, class Container>
bool join(const Polygon_with_holes_2<Kernel, Container> & p2,
          const Polygon_with_holes_2<Kernel, Container> & p1,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container> > & res);


/*!
 writes the union of the general polygons `p1` and `p2` into the polygon with holes `res`.
 Returns `true` if the two given polygons overlap.
 */
template <class Traits>
bool join(const General_polygon_2<Traits> & p1,
          const General_polygon_2<Traits> & p2,
          General_polygon_with_holes_2<General_polygon_2<Traits> > & res);


/*!
 writes the union of the polygons `p1` and `p2` into the polygon with holes `res`.
 Returns `true` if the two given polygons overlap.
 */
template <class Traits>
bool join(const General_polygon_2<Traits> & p1,
          const General_polygon_with_holes_2<General_polygon_2<Traits> > & p2,
          General_polygon_with_holes_2<General_polygon_2<Traits> > & res);


/*!
 writes the union of the general polygons `p1` and `p2` into the polygon with holes `res`.
 Returns `true` if the two given polygons overlap.
 */
template <class Traits>
bool join(const General_polygon_with_holes_2<General_polygon_2<Traits> > & p2,
          const General_polygon_2<Traits> & p1,
          General_polygon_with_holes_2<General_polygon_2<Traits> > & res);


/*!
 writes the union of the general polygons `p1` and `p2` into the polygon with holes `res`.
 Returns `true` if the two given polygons overlap.
 */
template <class Polygon>
bool join(const General_polygon_with_holes_2<Polygon> & p1,
          const General_polygon_with_holes_2<Polygon> & p2,
          Traits::Polygon_with_holes_2 & res);


/*!
  computes the union of the general polygons (or general polygons with 
  holes) in the given range. (The value type of the input iterator is 
  used to distinguish between the two.) The result, represented by a set 
  of general polygon with holes, is written into the output iterator `oi`. 
  The output iterator is 
  returned. The value type of the `OutputIterator` is 
  `Traits::Polygon_with_holes_2`. 
*/
template <class InputIterator, class OutputIterator>
OutputIterator join(InputIterator begin, InputIterator end,
OutputIterator oi);

  /*!
    computes the union of the general polygons and general polygons 
    with holes in the given two ranges. The result, represented by a set 
    of general polygon with holes, is written into the output iterator `oi`. 
    The output iterator is 
    returned. The value type of the `OutputIterator` is 
    `Traits::Polygon_with_holes_2`. 
  */
template <class InputIterator1, class InputIterator2,
class OutputIterator>
OutputIterator join(InputIterator1 pgn_begin1, InputIterator1 pgn_end1,
InputIterator2 pgn_begin2, InputIterator2 pgn_end2,
OutputIterator oi);

/// @}
} /* namespace CGAL */

namespace CGAL {
/*!
\addtogroup boolean_oriented_side Oriented Side Functions
\ingroup PkgBooleanSetOperations2
\anchor ref_bso_oriented_side

Each one of these functions returns `ON_POSITIVE_SIDE` if the two 
given polygons `p1` and `p2` intersect in their interior, 
`ON_NEGATIVE_SIDE` if `p1` and `p2` do not intersect at 
all, and `ON_ORIENTED_BOUNDARY` if `p1` and `p2` intersect 
only in their boundaries. 

The signature of the function is
\code
  Oriented_side oriented_side(const Type1 & p1, const Type2 & p2);
\endcode

\cgalHeading{Parameters}

The types of the paramters of the `oriented_side()` function are any of the following combinations.

<div align="left"> 
<table cellpadding=3 border="1"> 
<tr><th>Type1</th><th>Type2</th></tr> 
<tr><td valign="center">Polygon_2</td><td valign="center">Polygon_2</td></tr> 
<tr><td valign="center">Polygon_2</td><td valign="center">Polygon_with_holes_2</td></tr> 
<tr><td valign="center">Polygon_with_holes_2</td><td valign="center">Polygon_2</td></tr> 
<tr><td valign="center">Polygon_with_holes_2</td><td valign="center">Polygon_with_holes_2</td></tr> 
<tr><td valign="center">General_polygon_2</td><td valign="center">General_polygon_2</td></tr> 
<tr><td valign="center">General_polygon_2</td><td valign="center">General_polygon_with_holes_2</td></tr> 
<tr><td valign="center">General_polygon_with_holes_2</td><td valign="center">General_polygon_2</td></tr> 
<tr><td valign="center">General_polygon_with_holes_2</td><td valign="center">General_polygon_with_holes_2</td></tr> 
</table> 
</div> 

\sa \link boolean_do_intersect `CGAL::do_intersect()` \endlink

*/


/// @{

template <class Kernel, class Container>
Oriented_side oriented_side(const Polygon_2<Kernel, Container> & p1,
                            const Polygon_2<Kernel, Container> & p2);


template <class Kernel, class Container>
Oriented_side oriented_side(const Polygon_2<Kernel, Container> & p1,
                            const Polygon_with_holes_2<Kernel, Container> & p2);


template <class Kernel, class Container>
Oriented_side oriented_side(const Polygon_with_holes_2<Kernel, Container> & p1,
                            const Polygon_2<Kernel, Container> & p2);


template <class Kernel, class Container>
Oriented_side oriented_side(const Polygon_with_holes_2<Kernel, Container> & p1,
                            const Polygon_with_holes_2<Kernel, Container> & p2);


template <class Traits>
Oriented_side oriented_side(const General_polygon_2<Traits> & p1,
                            const General_polygon_2<Traits> & p2);


template <class Traits>
Oriented_side oriented_side(const General_polygon_2<Traits> & p1,
                            const General_polygon_with_holes_2<General_polygon_2<Traits> > & p2);


template <class Traits>
Oriented_side oriented_side(const General_polygon_with_holes_2<General_polygon_2<Traits> > & p1,
                            const General_polygon_2<Traits> & p2);


template <class Polygon>
Oriented_side oriented_side(const General_polygon_with_holes_2<Polygon> & p1,
                            const General_polygon_with_holes_2<Polygon> & p2);

/// @}
} /* namespace CGAL */

namespace CGAL {

/*!
\addtogroup boolean_symmetric_difference Symmetric Difference Functions
\ingroup PkgBooleanSetOperations2
\anchor ref_bso_symmetric_difference 

Each one of these functions computes the symmetric difference between 
two given polygons `p1` and `p2`, and inserts the resulting 
polygons with holes into an output container through the output 
iterator `oi`. The value type of the `OutputIterator` is either 
`Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

The signature of the function is
\code
  OutputIterator symmetric_difference(const Type1 & p1, const Type2 & p2, OutputIterator oi);
\endcode

\cgalHeading{Parameters}

The types of the paramters of the `symmetric_difference()` function are any of the following combinations.

<div align="left"> 
<table cellpadding=3 border="1"> 
<tr><th> Arg 1 type</th><th>Arg 2 type</th></tr> 
<tr><td valign="center">Polygon_2</td><td valign="center">Polygon_2</td></tr> 
<tr><td valign="center">Polygon_2</td><td valign="center">Polygon_with_holes_2</td></tr> 
<tr><td valign="center">Polygon_with_holes_2</td><td valign="center">Polygon_2</td></tr> 
<tr><td valign="center">Polygon_with_holes_2</td><td valign="center">Polygon_with_holes_2</td></tr> 
<tr><td valign="center">General_polygon_2</td><td valign="center">General_polygon_2</td></tr> 
<tr><td valign="center">General_polygon_2</td><td valign="center">General_polygon_with_holes_2</td></tr> 
<tr><td valign="center">General_polygon_with_holes_2</td><td valign="center">General_polygon_2</td></tr> 
<tr><td valign="center">General_polygon_with_holes_2</td><td valign="center">General_polygon_with_holes_2</td></tr> 
</table> 
</div> 


\sa \link boolean_do_intersect `CGAL::do_intersect()` \endlink
\sa \link boolean_intersection `CGAL::intersection()` \endlink
\sa \link boolean_join `CGAL::join()` \endlink
\sa \link boolean_difference `CGAL::difference()` \endlink

*/

/// @{

OutputIterator symmetric_difference(const Type1 & p1, const Type2 & p2,
                                    OutputIterator oi);


template <class Kernel, class Container, class OutputIterator>
OutputIterator symmetric_difference(const Polygon_2<Kernel, Container> & p1,
                                    const Polygon_2<Kernel, Container> & p2,
                                    OutputIterator oi);


template <class Kernel, class Container, class OutputIterator>
OutputIterator 
symmetric_difference(const Polygon_2<Kernel, Container> & p1,
                     const Polygon_with_holes_2<Kernel, Container> & p2,
                     OutputIterator oi);

template <class Kernel, class Container, class OutputIterator>
OutputIterator 
symmetric_difference(const Polygon_with_holes_2<Kernel, Container> & p1,
                     const Polygon_2<Kernel, Container> & p2,
                     OutputIterator oi);


template <class Kernel, class Container, class OutputIterator>
OutputIterator 
symmetric_difference(const Polygon_with_holes_2<Kernel, Container> & p1,
                     const Polygon_with_holes_2<Kernel, Container> & p2,
                     OutputIterator oi);


template <class Traits, class OutputIterator>
OutputIterator symmetric_difference(const General_polygon_2<Traits> & p1,
                                    const General_polygon_2<Traits> & p2,
                                    OutputIterator oi);


template <class Traits, class OutputIterator>
OutputIterator 
symmetric_difference(const General_polygon_with_holes_2<General_polygon_2<Traits> > & p1,
                     const General_polygon_2<Traits> & p2,
                     OutputIterator oi);


template <class Traits, class OutputIterator>
OutputIterator 
symmetric_difference(const General_polygon_2<Traits> & p1,
                     const General_polygon_with_holes_2<General_polygon_2<Traits> > & p2,
                     OutputIterator oi);


template <class Polygon, class OutputIterator>
OutputIterator 
symmetric_difference(const General_polygon_with_holes_2<Polygon> & p1,
                     const General_polygon_with_holes_2<Polygon> & p2,
                     OutputIterator oi);

  /*!
    computes the symmetric difference of the general polygons (or general 
    polygons with holes) in the given range. A point is contained in the 
    symmetric difference, if and only if it is contained in an odd number of 
    input polygons. (The value type of the input iterator is used to 
    distinguish between the two.) The result, represented by a set 
    of general polygon with holes, is inserted into an output container 
    through a given output iterator `oi`. The output iterator is 
    returned. The value type of the `OutputIterator` is 
    `Traits::Polygon_with_holes_2`. 
  */
template <class InputIterator, class OutputIterator>
OutputIterator symmetric_difference(InputIterator begin, InputIterator end,
                                    OutputIterator oi);

  /*!
    computes the symmetric difference of the general polygons and general polygons 
    with holes in the given two ranges. A point is contained in the 
    symmetric difference, if and only if it is contained in an odd number of 
    input polygons. The result, represented by a set of general polygon with 
    holes, is inserted into an output container through a given output 
    iterator `oi`. The output iterator is returned. The value type of 
    the `OutputIterator` is `Traits::Polygon_with_holes_2`. 
  */
template <class InputIterator1, class InputIterator2, class OutputIterator>
OutputIterator symmetric_difference(InputIterator1 pgn_begin1, 
                                    InputIterator1 pgn_end1,
                                    InputIterator2 pgn_begin2, 
                                    InputIterator2 pgn_end2,
                                    OutputIterator oi);
/// @}

} /* namespace CGAL */

