namespace CGAL {

/*!
\addtogroup boolean_complement Complement Functions
\ingroup PkgBooleanSetOperations2
\anchor ref_bso_complement 

Each one of these functions computes the complement of a given 
polygon `pgn`, and stores the resulting polygon with holes in 
`res`. 

<div align="center"> 
<table cellpadding=3 border="1"> 
<tr><th> Arg type</th></tr> 
<tr><td valign="center">Polygon_2</td></tr> 
<tr><td valign="center">General_polygon_2</td></tr> 
</table> 
</div> 

Each one of these functions computes the complement of a given 
polygon `pgn`, inserts the resulting polygons with 
holes into an output container through a given output iterator 
`oi`, and returns the output iterator. The value type of the 
`OutputIterator` is either `Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

<div align="center"> 
<table cellpadding=3 border="1"> 
<tr><th> Arg type</th></tr> 
<tr><td valign="center">Polygon_with_holes_2</td></tr> 
<tr><td valign="center">General_polygon_with_holes_2</td></tr> 
</table> 
</div> 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 
*/
/// @{
template <class Kernel, class Container>
void complement(const Polygon_2<Kernel, Container> & pgn, Polygon_with_holes_2<Kernel, Container> & res);
template <class Traits>
void complement(const General_polygon_2<Traits> & pgn, General_polygon_with_holes_2<Traits> & res);
template <class Traits, class OutputIterator>
OutputIterator complement(const Polygon_with_holes_2<Kernel, Container> & pgn, OutputIterator oi);
template <class Traits, class OutputIterator>
OutputIterator complement(const General_polygon_with_holes_2<General_polygon_2<Traits> > & pgn, OutputIterator oi);
/// @}

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the difference between two given 
polygons `p1` and `p2`, and inserts the resulting polygons 
with holes into an output container through the output iterator `oi`. 
The value type of the `OutputIterator` is either 
`Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::symmetric_difference` 

*/

OutputIterator difference(const Type1 & p1, const Type2 & p2,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the difference between two given 
polygons `p1` and `p2`, and inserts the resulting polygons 
with holes into an output container through the output iterator `oi`. 
The value type of the `OutputIterator` is either 
`Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::symmetric_difference` 

*/
template <class Kernel, class Container, class OutputIterator>
OutputIterator difference(const Polygon_2<Kernel, Container> & p1,
const Polygon_2<Kernel, Container> & p2,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the difference between two given 
polygons `p1` and `p2`, and inserts the resulting polygons 
with holes into an output container through the output iterator `oi`. 
The value type of the `OutputIterator` is either 
`Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::symmetric_difference` 

*/
template <class Kernel, class Container, class OutputIterator>
OutputIterator difference(const Polygon_2<Kernel, Container> & p1,
const Polygon_with_holes_2<Kernel,Container> & p2,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the difference between two given 
polygons `p1` and `p2`, and inserts the resulting polygons 
with holes into an output container through the output iterator `oi`. 
The value type of the `OutputIterator` is either 
`Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::symmetric_difference` 

*/
template <class Kernel, class Container, class OutputIterator>
OutputIterator difference(const Polygon_with_holes_2<Kernel, Container> & p1,
const Polygon_2<Kernel, Container> & p2,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the difference between two given 
polygons `p1` and `p2`, and inserts the resulting polygons 
with holes into an output container through the output iterator `oi`. 
The value type of the `OutputIterator` is either 
`Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::symmetric_difference` 

*/
template <class Kernel, class Container, class OutputIterator>
OutputIterator difference(const Polygon_with_holes_2<Kernel, Container> & p1,
const Polygon_with_holes_2<Kernel, Container> & p2,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the difference between two given 
polygons `p1` and `p2`, and inserts the resulting polygons 
with holes into an output container through the output iterator `oi`. 
The value type of the `OutputIterator` is either 
`Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::symmetric_difference` 

*/
template <class Traits, class OutputIterator>
OutputIterator difference(const General_polygon_2<Traits> & p1,
const General_polygon_2<Traits> & p2,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the difference between two given 
polygons `p1` and `p2`, and inserts the resulting polygons 
with holes into an output container through the output iterator `oi`. 
The value type of the `OutputIterator` is either 
`Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::symmetric_difference` 

*/
template <class Traits, class OutputIterator>
OutputIterator difference(const General_polygon_with_holes_2<General_polygon_2<Traits> > & p1,
const General_polygon_2<Traits> & p2,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the difference between two given 
polygons `p1` and `p2`, and inserts the resulting polygons 
with holes into an output container through the output iterator `oi`. 
The value type of the `OutputIterator` is either 
`Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::symmetric_difference` 

*/
template <class Traits, class OutputIterator>
OutputIterator difference(const General_polygon_2<Traits> & p1,
const General_polygon_with_holes_2<General_polygon_2<Traits> > & p2,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the difference between two given 
polygons `p1` and `p2`, and inserts the resulting polygons 
with holes into an output container through the output iterator `oi`. 
The value type of the `OutputIterator` is either 
`Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::symmetric_difference` 

*/
template <class Polygon, class OutputIterator>
OutputIterator difference(const General_polygon_with_holes_2<Polygon> & p1,
const General_polygon_with_holes_2<Polygon> & p2);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions return `true`, if the two given polygons 
`p1` and `p2` intersect in their interior, and `false` 
otherwise. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Returns `true`, if the set of general polygons (or general 
polygons with holes) in the given range intersect in their interior, 
and `false` otherwise. (The value type of the input iterator is 
used to distinguish between the two). 

Returns `true`, if the set of general polygons and general 
polygons with holes in the given two ranges respectively intersect in 
their interior, and `false` otherwise. 

\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
bool do_intersect(const Type1 & p1, const Type2 & p2);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions return `true`, if the two given polygons 
`p1` and `p2` intersect in their interior, and `false` 
otherwise. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Returns `true`, if the set of general polygons (or general 
polygons with holes) in the given range intersect in their interior, 
and `false` otherwise. (The value type of the input iterator is 
used to distinguish between the two). 

Returns `true`, if the set of general polygons and general 
polygons with holes in the given two ranges respectively intersect in 
their interior, and `false` otherwise. 

\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class Kernel, class Container>
bool do_intersect(const Polygon_2<Kernel, Container> & p1,
const Polygon_2<Kernel, Container> & p2);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions return `true`, if the two given polygons 
`p1` and `p2` intersect in their interior, and `false` 
otherwise. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Returns `true`, if the set of general polygons (or general 
polygons with holes) in the given range intersect in their interior, 
and `false` otherwise. (The value type of the input iterator is 
used to distinguish between the two). 

Returns `true`, if the set of general polygons and general 
polygons with holes in the given two ranges respectively intersect in 
their interior, and `false` otherwise. 

\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class Kernel, class Container>
bool do_intersect(const Polygon_2<Kernel, Container> & p1,
const Polygon_with_holes_2<Kernel, Container> & p2);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions return `true`, if the two given polygons 
`p1` and `p2` intersect in their interior, and `false` 
otherwise. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Returns `true`, if the set of general polygons (or general 
polygons with holes) in the given range intersect in their interior, 
and `false` otherwise. (The value type of the input iterator is 
used to distinguish between the two). 

Returns `true`, if the set of general polygons and general 
polygons with holes in the given two ranges respectively intersect in 
their interior, and `false` otherwise. 

\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class Kernel, class Container>
bool do_intersect(const Polygon_with_holes_2<Kernel, Container> & p1,
const Polygon_2<Kernel, Container> & p2);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions return `true`, if the two given polygons 
`p1` and `p2` intersect in their interior, and `false` 
otherwise. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Returns `true`, if the set of general polygons (or general 
polygons with holes) in the given range intersect in their interior, 
and `false` otherwise. (The value type of the input iterator is 
used to distinguish between the two). 

Returns `true`, if the set of general polygons and general 
polygons with holes in the given two ranges respectively intersect in 
their interior, and `false` otherwise. 

\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class Kernel, class Container>
bool do_intersect(const Polygon_with_holes_2<Kernel, Container> & p1,
const Polygon_with_holes_2<Kernel, Container> & p2);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions return `true`, if the two given polygons 
`p1` and `p2` intersect in their interior, and `false` 
otherwise. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Returns `true`, if the set of general polygons (or general 
polygons with holes) in the given range intersect in their interior, 
and `false` otherwise. (The value type of the input iterator is 
used to distinguish between the two). 

Returns `true`, if the set of general polygons and general 
polygons with holes in the given two ranges respectively intersect in 
their interior, and `false` otherwise. 

\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class Traits>
bool do_intersect(const General_polygon_2<Traits> & p1,
const General_polygon_2<Traits> & p2);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions return `true`, if the two given polygons 
`p1` and `p2` intersect in their interior, and `false` 
otherwise. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Returns `true`, if the set of general polygons (or general 
polygons with holes) in the given range intersect in their interior, 
and `false` otherwise. (The value type of the input iterator is 
used to distinguish between the two). 

Returns `true`, if the set of general polygons and general 
polygons with holes in the given two ranges respectively intersect in 
their interior, and `false` otherwise. 

\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class Traits>
bool do_intersect(const General_polygon_2<Traits> & p1,
const General_polygon_with_holes_2<General_polygon_2<Traits> > & p2);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions return `true`, if the two given polygons 
`p1` and `p2` intersect in their interior, and `false` 
otherwise. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Returns `true`, if the set of general polygons (or general 
polygons with holes) in the given range intersect in their interior, 
and `false` otherwise. (The value type of the input iterator is 
used to distinguish between the two). 

Returns `true`, if the set of general polygons and general 
polygons with holes in the given two ranges respectively intersect in 
their interior, and `false` otherwise. 

\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class Traits>
bool do_intersect(const General_polygon_with_holes_2<General_polygon_2<Traits> > & p1,
const General_polygon_2<Traits> & p2);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions return `true`, if the two given polygons 
`p1` and `p2` intersect in their interior, and `false` 
otherwise. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Returns `true`, if the set of general polygons (or general 
polygons with holes) in the given range intersect in their interior, 
and `false` otherwise. (The value type of the input iterator is 
used to distinguish between the two). 

Returns `true`, if the set of general polygons and general 
polygons with holes in the given two ranges respectively intersect in 
their interior, and `false` otherwise. 

\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class Polygon>
bool do_intersect(const General_polygon_with_holes_2<Polygon> & p1,
const General_polygon_with_holes_2<Polygon> & p2);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions return `true`, if the two given polygons 
`p1` and `p2` intersect in their interior, and `false` 
otherwise. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Returns `true`, if the set of general polygons (or general 
polygons with holes) in the given range intersect in their interior, 
and `false` otherwise. (The value type of the input iterator is 
used to distinguish between the two). 

Returns `true`, if the set of general polygons and general 
polygons with holes in the given two ranges respectively intersect in 
their interior, and `false` otherwise. 

\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class InputIterator>
bool do_intersect(InputIterator begin, InputIterator end);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions return `true`, if the two given polygons 
`p1` and `p2` intersect in their interior, and `false` 
otherwise. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Returns `true`, if the set of general polygons (or general 
polygons with holes) in the given range intersect in their interior, 
and `false` otherwise. (The value type of the input iterator is 
used to distinguish between the two). 

Returns `true`, if the set of general polygons and general 
polygons with holes in the given two ranges respectively intersect in 
their interior, and `false` otherwise. 

\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class InputIterator1, class InputIterator2>
bool do_intersect(InputIterator1 pgn_begin1,
InputIterator1 pgn_end1,
InputIterator2 pgn_begin2,
InputIterator2 pgn_end2);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the intersection of two given 
polygons `p1` and `p2`, inserts the resulting polygons with 
holes into an output container through a given output iterator 
`oi`, and returns the output iterator. The value type of the 
`OutputIterator` is either `Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Computes the intersection of the general polygons (or general polygons with 
holes) in the given range. (The value type of the input iterator is 
used to distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the intersection of the general polygons and general polygons 
with holes in the given two ranges. The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::join` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/

OutputIterator intersection(const Type1 & p1, const Type2 & p2,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the intersection of two given 
polygons `p1` and `p2`, inserts the resulting polygons with 
holes into an output container through a given output iterator 
`oi`, and returns the output iterator. The value type of the 
`OutputIterator` is either `Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Computes the intersection of the general polygons (or general polygons with 
holes) in the given range. (The value type of the input iterator is 
used to distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the intersection of the general polygons and general polygons 
with holes in the given two ranges. The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::join` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class Kernel, class Container, class OutputIterator>
OutputIterator intersection(const Polygon_2<Kernel, Container> & p1,
const Polygon_2<Kernel, Container> & p2,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the intersection of two given 
polygons `p1` and `p2`, inserts the resulting polygons with 
holes into an output container through a given output iterator 
`oi`, and returns the output iterator. The value type of the 
`OutputIterator` is either `Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Computes the intersection of the general polygons (or general polygons with 
holes) in the given range. (The value type of the input iterator is 
used to distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the intersection of the general polygons and general polygons 
with holes in the given two ranges. The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::join` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class Kernel, class Container, class OutputIterator>
OutputIterator intersection(const Polygon_2<Kernel, Container> & p1,
const Polygon_with_holes_2<Kernel, Container> & p2,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the intersection of two given 
polygons `p1` and `p2`, inserts the resulting polygons with 
holes into an output container through a given output iterator 
`oi`, and returns the output iterator. The value type of the 
`OutputIterator` is either `Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Computes the intersection of the general polygons (or general polygons with 
holes) in the given range. (The value type of the input iterator is 
used to distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the intersection of the general polygons and general polygons 
with holes in the given two ranges. The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::join` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class Kernel, class Container, class OutputIterator>
OutputIterator intersection(const Polygon_with_holes_2<Kernel, Container> & p1,
const Polygon_2<Kernel, Container> & p2,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the intersection of two given 
polygons `p1` and `p2`, inserts the resulting polygons with 
holes into an output container through a given output iterator 
`oi`, and returns the output iterator. The value type of the 
`OutputIterator` is either `Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Computes the intersection of the general polygons (or general polygons with 
holes) in the given range. (The value type of the input iterator is 
used to distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the intersection of the general polygons and general polygons 
with holes in the given two ranges. The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::join` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class Kernel, class Container, class OutputIterator>
OutputIterator intersection(const Polygon_with_holes_2<Kernel, Container> & p1,
const Polygon_with_holes_2<Kernel, Container> & p2,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the intersection of two given 
polygons `p1` and `p2`, inserts the resulting polygons with 
holes into an output container through a given output iterator 
`oi`, and returns the output iterator. The value type of the 
`OutputIterator` is either `Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Computes the intersection of the general polygons (or general polygons with 
holes) in the given range. (The value type of the input iterator is 
used to distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the intersection of the general polygons and general polygons 
with holes in the given two ranges. The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::join` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class Traits, class OutputIterator>
OutputIterator intersection(const General_polygon_2<Traits> & p1,
const General_polygon_2<Traits> & p2,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the intersection of two given 
polygons `p1` and `p2`, inserts the resulting polygons with 
holes into an output container through a given output iterator 
`oi`, and returns the output iterator. The value type of the 
`OutputIterator` is either `Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Computes the intersection of the general polygons (or general polygons with 
holes) in the given range. (The value type of the input iterator is 
used to distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the intersection of the general polygons and general polygons 
with holes in the given two ranges. The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::join` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class Traits, class OutputIterator>
OutputIterator intersection(const General_polygon_with_holes_2<General_polygon_2<Traits> > & p1,
const General_polygon_2<Traits> & p2,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the intersection of two given 
polygons `p1` and `p2`, inserts the resulting polygons with 
holes into an output container through a given output iterator 
`oi`, and returns the output iterator. The value type of the 
`OutputIterator` is either `Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Computes the intersection of the general polygons (or general polygons with 
holes) in the given range. (The value type of the input iterator is 
used to distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the intersection of the general polygons and general polygons 
with holes in the given two ranges. The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::join` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class Traits, class OutputIterator>
OutputIterator intersection(const General_polygon_2<Traits> & p1,
const General_polygon_with_holes_2<General_polygon_2<Traits> > & p2,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the intersection of two given 
polygons `p1` and `p2`, inserts the resulting polygons with 
holes into an output container through a given output iterator 
`oi`, and returns the output iterator. The value type of the 
`OutputIterator` is either `Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Computes the intersection of the general polygons (or general polygons with 
holes) in the given range. (The value type of the input iterator is 
used to distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the intersection of the general polygons and general polygons 
with holes in the given two ranges. The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::join` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class Polygon, class OutputIterator>
OutputIterator intersection(const General_polygon_with_holes_2<Polygon> & p1,
const General_polygon_with_holes_2<Polygon> & p2,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the intersection of two given 
polygons `p1` and `p2`, inserts the resulting polygons with 
holes into an output container through a given output iterator 
`oi`, and returns the output iterator. The value type of the 
`OutputIterator` is either `Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Computes the intersection of the general polygons (or general polygons with 
holes) in the given range. (The value type of the input iterator is 
used to distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the intersection of the general polygons and general polygons 
with holes in the given two ranges. The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::join` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class InputIterator, class OutputIterator>
OutputIterator intersection(InputIterator begin, InputIterator end,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the intersection of two given 
polygons `p1` and `p2`, inserts the resulting polygons with 
holes into an output container through a given output iterator 
`oi`, and returns the output iterator. The value type of the 
`OutputIterator` is either `Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Computes the intersection of the general polygons (or general polygons with 
holes) in the given range. (The value type of the input iterator is 
used to distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the intersection of the general polygons and general polygons 
with holes in the given two ranges. The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::join` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class InputIterator1, class InputIterator2,
class OutputIterator>
OutputIterator intersection(InputIterator1 pgn_begin1,
InputIterator1 pgn_end1,
InputIterator2 pgn_begin2,
InputIterator2 pgn_end2,
OutputIterator oi);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the union of two given polygons 
`p1` and `p2`. If the two given polygons overlap, it returns 
`true`, and places the resulting polygon in `p`. Otherwise, it 
returns `false`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
<table cellpadding=3 border="1"> 
<tr><th> Arg 1 type</th><th>Arg 2 type</th></tr> 
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
CONVERROR EndHtmlOnly 

Computes the union of the general polygons (or general polygons with 
holes) in the given range. (The value type of the input iterator is 
used to distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the union of the general polygons and general polygons 
with holes in the given two ranges. The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/

bool join(const Type1 & p1, const Type2 & p2,
General_polygon_with_holes_2 & p);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the union of two given polygons 
`p1` and `p2`. If the two given polygons overlap, it returns 
`true`, and places the resulting polygon in `p`. Otherwise, it 
returns `false`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
<table cellpadding=3 border="1"> 
<tr><th> Arg 1 type</th><th>Arg 2 type</th></tr> 
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
CONVERROR EndHtmlOnly 

Computes the union of the general polygons (or general polygons with 
holes) in the given range. (The value type of the input iterator is 
used to distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the union of the general polygons and general polygons 
with holes in the given two ranges. The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class Kernel, class Container>
bool join(const Polygon_2<Kernel, Container> & p1,
const Polygon_2<Kernel, Container> & p2,
General_polygon_with_holes_2<Polygon_2<Kernel, Container> > & p);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the union of two given polygons 
`p1` and `p2`. If the two given polygons overlap, it returns 
`true`, and places the resulting polygon in `p`. Otherwise, it 
returns `false`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
<table cellpadding=3 border="1"> 
<tr><th> Arg 1 type</th><th>Arg 2 type</th></tr> 
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
CONVERROR EndHtmlOnly 

Computes the union of the general polygons (or general polygons with 
holes) in the given range. (The value type of the input iterator is 
used to distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the union of the general polygons and general polygons 
with holes in the given two ranges. The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class Kernel, class Container>
bool join(const Polygon_2<Kernel, Container> & p1,
const Polygon_with_holes_2<Kernel,Container> & p2,
General_polygon_with_holes_2<Polygon_2<Kernel, Container> > & p);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the union of two given polygons 
`p1` and `p2`. If the two given polygons overlap, it returns 
`true`, and places the resulting polygon in `p`. Otherwise, it 
returns `false`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
<table cellpadding=3 border="1"> 
<tr><th> Arg 1 type</th><th>Arg 2 type</th></tr> 
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
CONVERROR EndHtmlOnly 

Computes the union of the general polygons (or general polygons with 
holes) in the given range. (The value type of the input iterator is 
used to distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the union of the general polygons and general polygons 
with holes in the given two ranges. The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class Kernel, class Container>
bool join(const Polygon_with_holes_2<Kernel, Container> & p2,
const Polygon_2<Kernel, Container> & p1,
General_polygon_with_holes_2<Polygon_2<Kernel, Container> > & p);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the union of two given polygons 
`p1` and `p2`. If the two given polygons overlap, it returns 
`true`, and places the resulting polygon in `p`. Otherwise, it 
returns `false`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
<table cellpadding=3 border="1"> 
<tr><th> Arg 1 type</th><th>Arg 2 type</th></tr> 
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
CONVERROR EndHtmlOnly 

Computes the union of the general polygons (or general polygons with 
holes) in the given range. (The value type of the input iterator is 
used to distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the union of the general polygons and general polygons 
with holes in the given two ranges. The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class Kernel, class Container>
bool join(const Polygon_with_holes_2<Kernel, Container> & p2,
const Polygon_with_holes_2<Kernel, Container> & p1,
General_polygon_with_holes_2<Polygon_2<Kernel, Container> > & p);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the union of two given polygons 
`p1` and `p2`. If the two given polygons overlap, it returns 
`true`, and places the resulting polygon in `p`. Otherwise, it 
returns `false`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
<table cellpadding=3 border="1"> 
<tr><th> Arg 1 type</th><th>Arg 2 type</th></tr> 
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
CONVERROR EndHtmlOnly 

Computes the union of the general polygons (or general polygons with 
holes) in the given range. (The value type of the input iterator is 
used to distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the union of the general polygons and general polygons 
with holes in the given two ranges. The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class Traits>
bool join(const General_polygon_2<Traits> & p1,
const General_polygon_2<Traits> & p2,
General_polygon_with_holes_2<General_polygon_2<Traits> > & p);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the union of two given polygons 
`p1` and `p2`. If the two given polygons overlap, it returns 
`true`, and places the resulting polygon in `p`. Otherwise, it 
returns `false`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
<table cellpadding=3 border="1"> 
<tr><th> Arg 1 type</th><th>Arg 2 type</th></tr> 
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
CONVERROR EndHtmlOnly 

Computes the union of the general polygons (or general polygons with 
holes) in the given range. (The value type of the input iterator is 
used to distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the union of the general polygons and general polygons 
with holes in the given two ranges. The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class Traits>
bool join(const General_polygon_2<Traits> & p1,
const General_polygon_with_holes_2<General_polygon_2<Traits> > & p2,
General_polygon_with_holes_2<General_polygon_2<Traits> > & p);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the union of two given polygons 
`p1` and `p2`. If the two given polygons overlap, it returns 
`true`, and places the resulting polygon in `p`. Otherwise, it 
returns `false`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
<table cellpadding=3 border="1"> 
<tr><th> Arg 1 type</th><th>Arg 2 type</th></tr> 
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
CONVERROR EndHtmlOnly 

Computes the union of the general polygons (or general polygons with 
holes) in the given range. (The value type of the input iterator is 
used to distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the union of the general polygons and general polygons 
with holes in the given two ranges. The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class Traits>
bool join(const General_polygon_with_holes_2<General_polygon_2<Traits> > & p2,
const General_polygon_2<Traits> & p1,
General_polygon_with_holes_2<General_polygon_2<Traits> > & p);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the union of two given polygons 
`p1` and `p2`. If the two given polygons overlap, it returns 
`true`, and places the resulting polygon in `p`. Otherwise, it 
returns `false`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
<table cellpadding=3 border="1"> 
<tr><th> Arg 1 type</th><th>Arg 2 type</th></tr> 
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
CONVERROR EndHtmlOnly 

Computes the union of the general polygons (or general polygons with 
holes) in the given range. (The value type of the input iterator is 
used to distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the union of the general polygons and general polygons 
with holes in the given two ranges. The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class Polygon>
bool join(const General_polygon_with_holes_2<Polygon> & p1,
const General_polygon_with_holes_2<Polygon> & p2,
Traits::Polygon_with_holes_2 & p);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the union of two given polygons 
`p1` and `p2`. If the two given polygons overlap, it returns 
`true`, and places the resulting polygon in `p`. Otherwise, it 
returns `false`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
<table cellpadding=3 border="1"> 
<tr><th> Arg 1 type</th><th>Arg 2 type</th></tr> 
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
CONVERROR EndHtmlOnly 

Computes the union of the general polygons (or general polygons with 
holes) in the given range. (The value type of the input iterator is 
used to distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the union of the general polygons and general polygons 
with holes in the given two ranges. The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class InputIterator, class OutputIterator>
OutputIterator join(InputIterator begin, InputIterator end,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the union of two given polygons 
`p1` and `p2`. If the two given polygons overlap, it returns 
`true`, and places the resulting polygon in `p`. Otherwise, it 
returns `false`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
<table cellpadding=3 border="1"> 
<tr><th> Arg 1 type</th><th>Arg 2 type</th></tr> 
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
CONVERROR EndHtmlOnly 

Computes the union of the general polygons (or general polygons with 
holes) in the given range. (The value type of the input iterator is 
used to distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the union of the general polygons and general polygons 
with holes in the given two ranges. The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::difference` 
\sa `CGAL::symmetric_difference` 

*/
template <class InputIterator1, class InputIterator2,
class OutputIterator>
OutputIterator join(InputIterator1 pgn_begin1, InputIterator1 pgn_end1,
InputIterator2 pgn_begin2, InputIterator2 pgn_end2,
OutputIterator oi);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions returns `ON_POSITIVE_SIDE` if the two 
given polygons `p1` and `p2` intersect in their interior, 
`ON_NEGATIVE_SIDE` if `p1` and `p2` do not intersect at 
all, and `ON_ORIENTED_BOUNDARY` if `p1` and `p2` intersect 
only in their boundaries. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

\sa `CGAL::do_intersect` 

*/
Oriented_side oriented_side(const Type1 & p1, const Type2 & p2);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions returns `ON_POSITIVE_SIDE` if the two 
given polygons `p1` and `p2` intersect in their interior, 
`ON_NEGATIVE_SIDE` if `p1` and `p2` do not intersect at 
all, and `ON_ORIENTED_BOUNDARY` if `p1` and `p2` intersect 
only in their boundaries. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

\sa `CGAL::do_intersect` 

*/
template <class Kernel, class Container>
Oriented_side oriented_side(const Polygon_2<Kernel, Container> & p1,
const Polygon_2<Kernel, Container> & p2);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions returns `ON_POSITIVE_SIDE` if the two 
given polygons `p1` and `p2` intersect in their interior, 
`ON_NEGATIVE_SIDE` if `p1` and `p2` do not intersect at 
all, and `ON_ORIENTED_BOUNDARY` if `p1` and `p2` intersect 
only in their boundaries. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

\sa `CGAL::do_intersect` 

*/
template <class Kernel, class Container>
Oriented_side oriented_side(const Polygon_2<Kernel, Container> & p1,
const Polygon_with_holes_2<Kernel, Container> & p2);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions returns `ON_POSITIVE_SIDE` if the two 
given polygons `p1` and `p2` intersect in their interior, 
`ON_NEGATIVE_SIDE` if `p1` and `p2` do not intersect at 
all, and `ON_ORIENTED_BOUNDARY` if `p1` and `p2` intersect 
only in their boundaries. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

\sa `CGAL::do_intersect` 

*/
template <class Kernel, class Container>
Oriented_side oriented_side(const Polygon_with_holes_2<Kernel, Container> & p1,
const Polygon_2<Kernel, Container> & p2);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions returns `ON_POSITIVE_SIDE` if the two 
given polygons `p1` and `p2` intersect in their interior, 
`ON_NEGATIVE_SIDE` if `p1` and `p2` do not intersect at 
all, and `ON_ORIENTED_BOUNDARY` if `p1` and `p2` intersect 
only in their boundaries. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

\sa `CGAL::do_intersect` 

*/
template <class Kernel, class Container>
Oriented_side oriented_side(const Polygon_with_holes_2<Kernel, Container> & p1,
const Polygon_with_holes_2<Kernel, Container> & p2);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions returns `ON_POSITIVE_SIDE` if the two 
given polygons `p1` and `p2` intersect in their interior, 
`ON_NEGATIVE_SIDE` if `p1` and `p2` do not intersect at 
all, and `ON_ORIENTED_BOUNDARY` if `p1` and `p2` intersect 
only in their boundaries. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

\sa `CGAL::do_intersect` 

*/
template <class Traits>
Oriented_side oriented_side(const General_polygon_2<Traits> & p1,
const General_polygon_2<Traits> & p2);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions returns `ON_POSITIVE_SIDE` if the two 
given polygons `p1` and `p2` intersect in their interior, 
`ON_NEGATIVE_SIDE` if `p1` and `p2` do not intersect at 
all, and `ON_ORIENTED_BOUNDARY` if `p1` and `p2` intersect 
only in their boundaries. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

\sa `CGAL::do_intersect` 

*/
template <class Traits>
Oriented_side oriented_side(const General_polygon_2<Traits> & p1,
const General_polygon_with_holes_2<General_polygon_2<Traits> > & p2);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions returns `ON_POSITIVE_SIDE` if the two 
given polygons `p1` and `p2` intersect in their interior, 
`ON_NEGATIVE_SIDE` if `p1` and `p2` do not intersect at 
all, and `ON_ORIENTED_BOUNDARY` if `p1` and `p2` intersect 
only in their boundaries. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

\sa `CGAL::do_intersect` 

*/
template <class Traits>
Oriented_side oriented_side(const General_polygon_with_holes_2<General_polygon_2<Traits> > & p1,
const General_polygon_2<Traits> & p2);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions returns `ON_POSITIVE_SIDE` if the two 
given polygons `p1` and `p2` intersect in their interior, 
`ON_NEGATIVE_SIDE` if `p1` and `p2` do not intersect at 
all, and `ON_ORIENTED_BOUNDARY` if `p1` and `p2` intersect 
only in their boundaries. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

\sa `CGAL::do_intersect` 

*/
template <class Polygon>
Oriented_side oriented_side(const General_polygon_with_holes_2<Polygon> & p1,
const General_polygon_with_holes_2<Polygon> & p2);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the symmetric difference between 
two given polygons `p1` and `p2`, and inserts the resulting 
polygons with holes into an output container through the output 
iterator `oi`. The value type of the `OutputIterator` is either 
`Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Computes the symmetric difference of the general polygons (or general 
polygons with holes) in the given range. A point is contained in the 
symmetric difference, if and only if it is contained in odd number of 
input polygons. (The value type of the input iterator is used to 
distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the union of the general polygons and general polygons 
with holes in the given two ranges. A point is contained in the 
symmetric difference, if and only if it is contained in odd number of 
input polygons. The result, represented by a set of general polygon with 
holes, is inserted into an output container through a given output 
iterator `oi`. The output iterator is returned. The value type of 
the `OutputIterator` is `Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::difference` 

*/

OutputIterator intersection(const Type1 & p1, const Type2 & p2,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the symmetric difference between 
two given polygons `p1` and `p2`, and inserts the resulting 
polygons with holes into an output container through the output 
iterator `oi`. The value type of the `OutputIterator` is either 
`Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Computes the symmetric difference of the general polygons (or general 
polygons with holes) in the given range. A point is contained in the 
symmetric difference, if and only if it is contained in odd number of 
input polygons. (The value type of the input iterator is used to 
distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the union of the general polygons and general polygons 
with holes in the given two ranges. A point is contained in the 
symmetric difference, if and only if it is contained in odd number of 
input polygons. The result, represented by a set of general polygon with 
holes, is inserted into an output container through a given output 
iterator `oi`. The output iterator is returned. The value type of 
the `OutputIterator` is `Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::difference` 

*/
template <class Kernel, class Container, class OutputIterator>
OutputIterator symmetric_difference(const Polygon_2<Kernel, Container> & p1,
const Polygon_2<Kernel, Container> & p2,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the symmetric difference between 
two given polygons `p1` and `p2`, and inserts the resulting 
polygons with holes into an output container through the output 
iterator `oi`. The value type of the `OutputIterator` is either 
`Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Computes the symmetric difference of the general polygons (or general 
polygons with holes) in the given range. A point is contained in the 
symmetric difference, if and only if it is contained in odd number of 
input polygons. (The value type of the input iterator is used to 
distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the union of the general polygons and general polygons 
with holes in the given two ranges. A point is contained in the 
symmetric difference, if and only if it is contained in odd number of 
input polygons. The result, represented by a set of general polygon with 
holes, is inserted into an output container through a given output 
iterator `oi`. The output iterator is returned. The value type of 
the `OutputIterator` is `Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::difference` 

*/
template <class Kernel, class Container, class OutputIterator>
OutputIterator 
symmetric_difference(const Polygon_2<Kernel, Container> & p1,
const Polygon_with_holes_2<Kernel, Container> & p2,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the symmetric difference between 
two given polygons `p1` and `p2`, and inserts the resulting 
polygons with holes into an output container through the output 
iterator `oi`. The value type of the `OutputIterator` is either 
`Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Computes the symmetric difference of the general polygons (or general 
polygons with holes) in the given range. A point is contained in the 
symmetric difference, if and only if it is contained in odd number of 
input polygons. (The value type of the input iterator is used to 
distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the union of the general polygons and general polygons 
with holes in the given two ranges. A point is contained in the 
symmetric difference, if and only if it is contained in odd number of 
input polygons. The result, represented by a set of general polygon with 
holes, is inserted into an output container through a given output 
iterator `oi`. The output iterator is returned. The value type of 
the `OutputIterator` is `Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::difference` 

*/
template <class Kernel, class Container, class OutputIterator>
OutputIterator 
symmetric_difference(const Polygon_with_holes_2<Kernel, Container> & p1,
const Polygon_2<Kernel, Container> & p2,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the symmetric difference between 
two given polygons `p1` and `p2`, and inserts the resulting 
polygons with holes into an output container through the output 
iterator `oi`. The value type of the `OutputIterator` is either 
`Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Computes the symmetric difference of the general polygons (or general 
polygons with holes) in the given range. A point is contained in the 
symmetric difference, if and only if it is contained in odd number of 
input polygons. (The value type of the input iterator is used to 
distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the union of the general polygons and general polygons 
with holes in the given two ranges. A point is contained in the 
symmetric difference, if and only if it is contained in odd number of 
input polygons. The result, represented by a set of general polygon with 
holes, is inserted into an output container through a given output 
iterator `oi`. The output iterator is returned. The value type of 
the `OutputIterator` is `Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::difference` 

*/
template <class Kernel, class Container, class OutputIterator>
OutputIterator 
symmetric_difference(const Polygon_with_holes_2<Kernel, Container> & p1,
const Polygon_with_holes_2<Kernel, Container> & p2,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the symmetric difference between 
two given polygons `p1` and `p2`, and inserts the resulting 
polygons with holes into an output container through the output 
iterator `oi`. The value type of the `OutputIterator` is either 
`Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Computes the symmetric difference of the general polygons (or general 
polygons with holes) in the given range. A point is contained in the 
symmetric difference, if and only if it is contained in odd number of 
input polygons. (The value type of the input iterator is used to 
distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the union of the general polygons and general polygons 
with holes in the given two ranges. A point is contained in the 
symmetric difference, if and only if it is contained in odd number of 
input polygons. The result, represented by a set of general polygon with 
holes, is inserted into an output container through a given output 
iterator `oi`. The output iterator is returned. The value type of 
the `OutputIterator` is `Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::difference` 

*/
template <class Traits, class OutputIterator>
OutputIterator symmetric_difference(const General_polygon_2<Traits> & p1,
const General_polygon_2<Traits> & p2,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the symmetric difference between 
two given polygons `p1` and `p2`, and inserts the resulting 
polygons with holes into an output container through the output 
iterator `oi`. The value type of the `OutputIterator` is either 
`Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Computes the symmetric difference of the general polygons (or general 
polygons with holes) in the given range. A point is contained in the 
symmetric difference, if and only if it is contained in odd number of 
input polygons. (The value type of the input iterator is used to 
distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the union of the general polygons and general polygons 
with holes in the given two ranges. A point is contained in the 
symmetric difference, if and only if it is contained in odd number of 
input polygons. The result, represented by a set of general polygon with 
holes, is inserted into an output container through a given output 
iterator `oi`. The output iterator is returned. The value type of 
the `OutputIterator` is `Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::difference` 

*/
template <class Traits, class OutputIterator>
OutputIterator 
symmetric_difference(const General_polygon_with_holes_2<General_polygon_2<Traits> > & p1,
const General_polygon_2<Traits> & p2,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the symmetric difference between 
two given polygons `p1` and `p2`, and inserts the resulting 
polygons with holes into an output container through the output 
iterator `oi`. The value type of the `OutputIterator` is either 
`Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Computes the symmetric difference of the general polygons (or general 
polygons with holes) in the given range. A point is contained in the 
symmetric difference, if and only if it is contained in odd number of 
input polygons. (The value type of the input iterator is used to 
distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the union of the general polygons and general polygons 
with holes in the given two ranges. A point is contained in the 
symmetric difference, if and only if it is contained in odd number of 
input polygons. The result, represented by a set of general polygon with 
holes, is inserted into an output container through a given output 
iterator `oi`. The output iterator is returned. The value type of 
the `OutputIterator` is `Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::difference` 

*/
template <class Traits, class OutputIterator>
OutputIterator 
symmetric_difference(const General_polygon_2<Traits> & p1,
const General_polygon_with_holes_2<General_polygon_2<Traits> > & p2,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the symmetric difference between 
two given polygons `p1` and `p2`, and inserts the resulting 
polygons with holes into an output container through the output 
iterator `oi`. The value type of the `OutputIterator` is either 
`Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Computes the symmetric difference of the general polygons (or general 
polygons with holes) in the given range. A point is contained in the 
symmetric difference, if and only if it is contained in odd number of 
input polygons. (The value type of the input iterator is used to 
distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the union of the general polygons and general polygons 
with holes in the given two ranges. A point is contained in the 
symmetric difference, if and only if it is contained in odd number of 
input polygons. The result, represented by a set of general polygon with 
holes, is inserted into an output container through a given output 
iterator `oi`. The output iterator is returned. The value type of 
the `OutputIterator` is `Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::difference` 

*/
template <class Polygon, class OutputIterator>
OutputIterator 
symmetric_difference(const General_polygon_with_holes_2<Polygon> & p1,
const General_polygon_with_holes_2<Polygon> & p2,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the symmetric difference between 
two given polygons `p1` and `p2`, and inserts the resulting 
polygons with holes into an output container through the output 
iterator `oi`. The value type of the `OutputIterator` is either 
`Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Computes the symmetric difference of the general polygons (or general 
polygons with holes) in the given range. A point is contained in the 
symmetric difference, if and only if it is contained in odd number of 
input polygons. (The value type of the input iterator is used to 
distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the union of the general polygons and general polygons 
with holes in the given two ranges. A point is contained in the 
symmetric difference, if and only if it is contained in odd number of 
input polygons. The result, represented by a set of general polygon with 
holes, is inserted into an output container through a given output 
iterator `oi`. The output iterator is returned. The value type of 
the `OutputIterator` is `Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::difference` 

*/
template <class InputIterator, class OutputIterator>
OutputIterator symmetric_difference(InputIterator begin, InputIterator end,
OutputIterator oi);

/*!
\ingroup PkgBooleanSetOperations2

Each one of these functions computes the symmetric difference between 
two given polygons `p1` and `p2`, and inserts the resulting 
polygons with holes into an output container through the output 
iterator `oi`. The value type of the `OutputIterator` is either 
`Polygon_with_holes_2` or 
`General_polygon_with_holes_2`. 

CONVERROR HtmlOnly needs treatment 
<div align="center"> 
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
CONVERROR EndHtmlOnly 

Computes the symmetric difference of the general polygons (or general 
polygons with holes) in the given range. A point is contained in the 
symmetric difference, if and only if it is contained in odd number of 
input polygons. (The value type of the input iterator is used to 
distinguish between the two.) The result, represented by a set 
of general polygon with holes, is inserted into an output container 
through a given output iterator `oi`. The output iterator is 
returned. The value type of the `OutputIterator` is 
`Traits::Polygon_with_holes_2`. 

Computes the union of the general polygons and general polygons 
with holes in the given two ranges. A point is contained in the 
symmetric difference, if and only if it is contained in odd number of 
input polygons. The result, represented by a set of general polygon with 
holes, is inserted into an output container through a given output 
iterator `oi`. The output iterator is returned. The value type of 
the `OutputIterator` is `Traits::Polygon_with_holes_2`. 

\sa `CGAL::do_intersect` 
\sa `CGAL::intersection` 
\sa `CGAL::join` 
\sa `CGAL::difference` 

*/
template <class InputIterator1, class InputIterator2,
class OutputIterator>
OutputIterator symmetric_difference(InputIterator1 pgn_begin1, 
InputIterator1 pgn_end1,
InputIterator2 pgn_begin2, 
InputIterator2 pgn_end2,
OutputIterator oi);

} /* namespace CGAL */

