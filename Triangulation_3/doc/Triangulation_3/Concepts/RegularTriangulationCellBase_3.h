
/*!
\ingroup PkgTriangulation3Concepts
\cgalConcept

The regular triangulation of a set of weighted points does not 
necessarily 
have one vertex for each of the input points. Some of the input 
weighted points have no cell in the dual power diagrams 
and therefore do not correspond to a vertex of the regular 
triangulation. 
Those weighted points are said to be <I>hidden</I> points. 
A point which is hidden at a given time may appear later as a vertex of 
the regular triangulation upon removal on some other weighted point. 
Therefore, hidden points have to be stored somewhere. 
The regular triangulation stores those hidden points 
in its cells. 

A hidden point can appear as vertex of the triangulation 
only when the 
three dimensional cell where its point component is located 
(the cell which hides it) 
is removed. Therefore we decided to store 
in each cell of a regular triangulation 
the list of hidden points that are located in the face. 
Thus points hidden by a face are easily reinserted in the triangulation 
when the face is removed. 

The base cell of a regular triangulation 
has to be a model 
of the concept `RegularTriangulationCellBase_3`, which refines 
the concept `TriangulationCellBase_3` by adding 
in the cell a container to store hidden points 
and an operator to compute its weighted circumcenter. 

\cgalRefines `TriangulationCellBase_3`

\cgalHasModel CGAL::Regular_triangulation_cell_base_3 
\cgalHasModel CGAL::Regular_triangulation_cell_base_with_weighted_circumcenter_3

\sa `RegularTriangulationTraits_3`

*/

class RegularTriangulationCellBase_3 {
public:

/// \name Types 
/// @{

/*!
Must be the same as the point type `RegularTriangulationTraits_3::Weighted_point_3` 
defined by the geometric traits class of the triangulation. 
*/ 
typedef unspecified_type Weighted_point; 

/*!
Iterator of value type Point 
*/ 
typedef unspecified_type Point_iterator; 


/// @} 

/// \name Access Functions 
/// @{

/*!
Returns an iterator pointing to the first hidden point. 
*/ 
Point_iterator hidden_points_begin(); 

/*!
Returns a past-the-end iterator. 
*/ 
Point_iterator hidden_points_end(); 

/// @} 

/// \name Setting 
/// @{

/*!
Adds `p` to the set of hidden points of the cell. 
*/ 
void hide_point(const Point & p); 

/// @}

/// \name Access functions
/// @{
/*!
Returns the weighted circumcenter of the cell, with no weight. 
`RegularTriangulationTraits_3` is the geometric traits class of the triangulation.
This operator is required only when the dual functions are called.
*/ 
const RegularTriangulationTraits_3::Bare_point& weighted_circumcenter( 
const RegularTriangulationTraits_3&gt = RegularTriangulationTraits_3()) const; 
/// @} 


}; /* end RegularTriangulationCellBase_3 */

