
namespace CGAL {

/*!
\ingroup PkgMesh2


This class implements a 2D mesh generator. 


\tparam CDT must be a 2D constrained Delaunay triangulation, and type `CDT::Face` 
should be a model of the concept `DelaunayMeshFaceBase_2`.
The geometric traits class of the instance of `CDT` has to be 
a model of the concept `DelaunayMeshTraits_2`. 

\tparam Criteria must be a model of the concept 
`MeshingCriteria_2`. This traits class defines the shape and size 
criteria for the triangles of the mesh. `Criteria::Face_handle` has to 
be the same as `CDT::Face_handle`. 

\cgalHeading{Using This Class}

The constructor of the class `Delaunay_mesher_2` takes a reference to a `CDT` 
as an argument. A call to the refinement method `refine_mesh()` will 
refine the constrained Delaunay triangulation into a mesh satisfying the 
size and shape criteria specified in the traits class. Note that if, during 
the life time of the `Delaunay_mesher_2` object, the triangulation is externally 
modified, any further call to its member methods may crash. Consider 
constructing a new `Delaunay_mesher_2` object if the triangulation has been 
modified. 

\cgalHeading{Meshing Domain}

The domain to be mesh is defined by the constrained edges and a set of seed 
points. The constrained edges divides the plane into several connected 
components. The mesh domain is either the union of the bounded connected 
components including at least one seed, or the union of the bounded 
connected components that do no contain any seed. Note that the unbounded 
component of the plane is never meshed. 

\sa `CGAL::lloyd_optimize_mesh_2`

*/
template< typename CDT, typename Criteria >
class Delaunay_mesher_2 {
public:

/// \name Types 
/// @{

/*!
the geometric traits class. 
*/ 
typedef CDT::Geom_traits Geom_traits; 






/*!
const iterator over defined seeds. Its 
value type is `Geom_traits::Point_2`. 
*/ 
typedef unspecified_type Seeds_iterator; 

/// @} 


/// \name Creation 
/// @{

/*!
Create a new mesher, working on `t`, with meshing criteria 
`criteria`. 
*/ 
Delaunay_mesher_2(CDT& t, Criteria criteria = Criteria()); 





/// @} 


/// \name Seeds functions 
/// The following functions are used to define seeds. 
/// @{

/*!
Sets seeds to the empty set. All 
finite connected components of the constrained triangulation will be 
refined. 
*/ 
void clear_seeds (); 





/*!
Sets seeds to the sequence `[begin, end)`. If `mark==true`, the mesh domain 
is the union of the bounded connected 
components including at least one seed. If 
`mark==false`, the domain is the union of 
the bounded components including no seed. Note 
that the unbounded component of the plane is 
never meshed. 

\tparam InputIterator must be an input iterator with  value type `Geom_traits::Point_2`. 
*/ 
template<class InputIterator> 
void set_seeds(InputIterator begin, InputIterator end, 
const bool mark=false); 





/*!
Start of the seeds sequence. 
*/ 
Seeds_const_iterator seeds_begin () const; 





/*!
Past the end of the seeds sequence. 
*/ 
Seeds_const_iterator seeds_end () const; 





/// @} 

/*!
\name Meshing methods 
The function `set_criteria()` scans all faces to recalculate the list of 
<I>bad faces</I>, that are faces not conforming to the meshing criteria. 
This function actually has an optional argument that permits to prevent 
this recalculation. The filling of the list of bad faces can then be done 
by a call to `set_bad_faces`. 
*/
/// @{

/*!
Refines the constrained Delaunay triangulation into a mesh 
satisfying the criteria defined by the traits. 

*/ 
void refine_mesh(); 





/*!
Returns a const reference to the criteria traits object. 
*/ 
const Criteria& get_criteria(); 





/*!
Assigns `criteria` to the criteria traits object. 
*/ 
void set_criteria(Criteria criteria); 





/*!
Assigns `criteria` to the criteria traits object. If 
`recalculate_bad_faces` is `false`, the list of bad faces is 
let empty and the function `set_bad_faces()` should be called before 
`refine_mesh`. 
*/ 
void set_criteria(Criteria criteria, bool 
recalculate_bad_faces); 





/*!
This method permits to set the list of bad triangles 
directly, from the sequence `[begin, end)`, so that the 
algorithm will not scan the whole set of triangles to 
find bad ones. To use if there is a non-naive way to 
find bad triangles. 

\tparam InputIterator must be an input iterator with value type `Face_handle`. 
*/ 
template <class InputIterator> 
void set_bad_faces(InputIterator begin, 
InputIterator end); 





/// @} 


/// \name Step by Step Operations 
/// The `Delaunay_mesher_2` class allows, for debugging or demos, to play the 
/// meshing algorithm step by step, using the following methods. 
/// @{

/*!
This method must be called just before the first 
call to the following step by step refinement method, 
that is when all vertices and constrained edges have been 
inserted into the constrained Delaunay triangulation. It 
must be called again before any subsequent calls to the 
step by step refinement method if new vertices or constrained 
edges have been inserted since the last call. 
*/ 
void init(); 





/*!
Tests if the step by step refinement algorithm is done. If it returns 
`true`, the following calls to `step_by_step_refine_mesh` will 
not insert any points, until some new constrained segments or points are 
inserted in the triangulation and `init` is called again. 
*/ 
bool is_refinement_done(); 





/*!
Applies one step of the algorithm, by inserting one point, if the 
algorithm is not done. Returns `false` iff no point has been inserted 
because the algorithm is done. 
*/ 
bool step_by_step_refine_mesh(); 





/// @}

}; /* end Delaunay_mesher_2 */

/*!
\ingroup PkgMesh2Functions

Refines the default domain defined by a constrained Delaunay
triangulation without seeds into a mesh satisfying the criteria
defined by the traits `criteria`. The domain of the mesh
covers all the connected components of the plane defined by the
constrained edges of `t`, except for the unbounded component.

\tparam CDT  must be 2D constrained Delaunay triangulation
and its geometric traits class must be a model of `DelaunayMeshTraits_2`.
The face of the constrained Delaunay triangulation must be a model of the concept `DelaunayMeshFaceBase_2`.

\tparam Criteria must be a model of the concept `MeshingCriteria_2`. 
`Criteria::Face_handle` must be the same as `CDT::Face_handle`.
*/
template<class CDT, class Criteria> 
void refine_Delaunay_mesh_2 (CDT &t, const Criteria& criteria = Criteria()); 


/*!

\ingroup PkgMesh2Functions

Refines the default domain defined by a constrained
Delaunay triangulation into a mesh
satisfying the criteria defined by the traits
`criteria`.The sequence `[begin, end)`
gives a set of seeds points, that defines the domain
to be meshed as follows. The constrained edges of
`t` partition the plane into connected components.
If `mark==true`, the mesh domain is the union of
the bounded connected components including at least
one seed. If `mark==false`, the domain is the
union of the bounded components including no seed.
Note that the unbounded component of the plane is
never meshed.

\tparam CDT  must be 2D constrained Delaunay triangulation
and its geometric traits class must be a model of `DelaunayMeshTraits_2`.
The face of the constrained Delaunay triangulation must be a model of the concept `DelaunayMeshFaceBase_2`.

\tparam Criteria must be a model of the concept `MeshingCriteria_2`. 
`Criteria::Face_handle` must be the same as `CDT::Face_handle`.
\tparam InputIterator must be an input iterator with value type `CDT::Geom_traits::Point_2`.


*/
template <class CDT, class Criteria, class InputIterator>
void refine_Delaunay_mesh_2(CDT& t,
InputIterator begin, InputIterator end,
const Criteria& criteria = Criteria(),
bool mark = false); 



} /* end namespace CGAL */
