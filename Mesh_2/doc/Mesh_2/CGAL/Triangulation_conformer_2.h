namespace CGAL {

/*!
\ingroup PkgMesh2Functions



Refines the constrained Delaunay triangulation `t` into a
conforming Delaunay triangulation. After a call to this function,
all edges of `t` are Delaunay edges. 

\tparam CDT must be a 2D constrained Delaunay triangulation
and its geometric traits class must be a model of `ConformingDelaunayTriangulationTraits_2`.
*/
template<class CDT> void make_conforming_Delaunay_2 (CDT &t); 

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgMesh2Functions



Refines the constrained Delaunay triangulation `t` into a
conforming Gabriel triangulation. After a call to this function, all
constrained edges of `t` have the <I>Gabriel property</I>: the
circle that has \f$ e\f$ as diameter does not contain any vertex from
the triangulation. 

\tparam CDT must be a 2D constrained Delaunay triangulation
and its geometric traits class must be a model of `ConformingDelaunayTriangulationTraits_2`.
*/
template<class CDT> void make_conforming_Gabriel_2 (CDT &t); 

} /* namespace CGAL */


namespace CGAL {

/*!
\ingroup PkgMesh2


The class `Triangulation_conformer_2` is an auxiliary class of 
`Delaunay_mesher_2<CDT>`. It permits to refine a constrained 
Delaunay triangulation into a conforming Delaunay or conforming 
Gabriel triangulation. For standard needs, consider using the global 
functions `make_conforming_Gabriel_2()` and 
`make_conforming_Delaunay_2()`. 


\tparam CDT must be a 2D constrained Delaunay triangulation
and its geometric traits class must be 
a model of the concept `ConformingDelaunayTriangulationTraits_2`. 

\cgalHeading{Using This Class}

The constructor of the class `Triangulation_conformer_2` takes a reference to a `CDT` 
as an argument. A call to the method `make_conforming_Delaunay()` or 
`make_conforming_Gabriel()` will refine this constrained Delaunay 
triangulation into a conforming Delaunay or conforming Gabriel 
triangulation. Note that if, during the life time of the `Triangulation_conformer_2` object, the triangulation is externally modified, any further call to its 
member methods may lead to undefined behavior. Consider reconstructing a 
new `Triangulation_conformer_2` object if the triangulation has been modified. 

The conforming methods insert points into constrained edges, thereby splitting 
them into several sub-constraints. You have access to the initial inserted 
constraints if you instantiate the template parameter by a 
`Constrained_triangulation_plus_2<CDT>`. 

*/
template< typename CDT >
class Triangulation_conformer_2 {
public:

/// \name Creation 
/// @{

/*!
Create a new conforming maker, working on `t`. 
*/ 
Triangulation_conformer_2(CDT& t); 





/// @} 


/// \name Conforming methods 
/// @{

/*!
Refines the triangulation into a conforming Delaunay triangulation. 
After a call to this method, all triangles fulfill the Delaunay property, 
that is the empty circle 
property. 
*/ 
void make_conforming_Delaunay(); 





/*!
Refines the triangulation into a conforming Gabriel triangulation. 
After a call to this method, all constrained edges \f$ e\f$ have the 
<I>Gabriel property</I>: the circle with diameter \f$ e\f$ 
does not contain any vertex of the triangulation. 
*/ 
void make_conforming_Gabriel(); 





/// @} 


/*!
\name Checking 
The following methods verify that the constrained triangulation is 
conforming Delaunay or conforming Gabriel. These methods scan the 
whole triangulation and their complexity is proportional to the number 
of edges. 
*/
/// @{

/*!
Returns `true` iff all triangles fulfill the Delaunay property. 
*/ 
bool is_conforming_Delaunay(); 





/*!
Returns `true` iff all constrained edges have the Gabriel property: 
their circumsphere is empty. 
*/ 
bool is_conforming_Gabriel(); 





/// @} 


/*!
\name Step by Step Operations 
The `Triangulation_conformer_2` class allows, for debugging or demos, to play the 
conforming algorithm step by step, using the following methods. They exist 
in two versions, depending on whether you want the triangulation to be 
conforming Delaunay or conforming Gabriel, respectively. Any call to a 
`step_by_step_conforming_XX` function requires a previous call to the 
corresponding function `init_XX` and Gabriel and Delaunay methods can 
not be mixed between two calls of `init_XX`. 
*/
/// @{

/*!
The method must be called after all points and constrained segments 
are inserted and before any call to the following methods. If some 
points or segments are then inserted in the triangulation, this 
method must be called again. 
*/ 
void init_Delaunay(); 





/*!
Applies one step of the algorithm, by inserting one point, if the 
algorithm is not done. Returns `false` iff no point has been inserted 
because the algorithm is done. 
*/ 
bool step_by_step_conforming_Delaunay (); 





/*!
Analog to 
`init_Delaunay` for Gabriel conforming. 
*/ 
void init_Gabriel(); 





/*!
Analog to 
`step_by_step_conforming_Delaunay()` for Gabriel conforming. 
*/ 
bool step_by_step_conforming_Gabriel (); 





/*!
Tests if the step by step conforming algorithm is done. If it 
returns `true`, the following calls to 
`step_by_step_conforming_XX` will not insert any points, until some 
new constrained segments or points are inserted in the triangulation and 
`init_XX` is called again. 
*/ 
bool is_conforming_done(); 





/// @}

}; /* end Triangulation_conformer_2 */
} /* end namespace CGAL */
