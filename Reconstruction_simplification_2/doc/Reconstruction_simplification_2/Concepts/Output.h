
/*!
\ingroup PkgReconstructionSimplification2Concepts
\cgalConcept

The OutputModel  refinement process involved in the 
function template `CGAL::make_surface_mesh()` 
is guided by a set of refinement criteria. 
The concept `SurfaceMeshFacetsCriteria_3` describes the type which 
handles those criteria. 
It corresponds to the requirements for the template parameter 
`FacetsCriteria` of the surface mesher function 
`CGAL::make_surface_mesh<SurfaceMeshC2T3,Surface,FacetsCriteria,Tag>()`. 

Typically the meshing criteria are a set 
of elementary criteria, each of which 
has to be met by the facets of the final mesh. 
The meshing algorithm eliminates in turn <I>bad</I> facets, i.e., 
facets that do not meet all the criteria. 

The size and quality of the final mesh 
depends on the order according to which bad facets 
are handled. Therefore, the meshing algorithm 
needs to be able to quantify the facet qualities and to compare 
the qualities of different faces. 
The type `SurfaceMeshFacetsCriteria_3::Quality` measures 
the quality of a mesh facet. 
Typically this quality 
is a multicomponent variable. Each component corresponds to 
one criterion and measures how much the facet deviates from 
meeting this criterion. Then, the comparison operator on qualities 
is just a lexicographical comparison. The meshing algorithm handles facets 
with lowest quality first. The qualities are computed by a function 
`is_bad(const Facet& f, const Quality& q)`. 

\cgalHasModel `CGAL::Surface_mesh_default_criteria_3<Tr>` 

\sa `CGAL::make_surface_mesh()` 

*/

class OutputModel {
public:

/// \name Types 
/// @{

/*!
The type of facets. This type has to match 
the `Facet` type in the triangulation type used by 
the mesher function. (This triangulation type 
is the type `SurfaceMeshC2T3::Triangulation` 
provided by the model of 
`SurfaceMeshComplex_2InTriangulation_3` plugged 
as first template parameter of 
`CGAL::make_surface_mesh()`). 
*/ 
typedef unspecified_type Facet; 

/*!
Default constructible, copy constructible, 
assignable, and less-than comparable type. 
*/ 
typedef unspecified_type Quality; 

/// @} 

/// \name Operations 
/// @{

/*!
Assigns the quality 
of the facet `f` to `q`, and returns `true` is `f` does 
not meet the criteria. 
*/ 
bool is_bad (const Facet& f, const Quality& q); 

/// @}

}; /* end SurfaceMeshFacetsCriteria_3 */

