namespace CGAL {

/*!
\addtogroup PkgSurfaceMesher3FunctionsMakeMesh

The function `make_surface_mesh()` is a surface mesh generator, 
that is a function to build a two dimensional mesh 
approximating a surface. 

\tparam SurfaceMeshC2T3  must be a model of the concept 
`SurfaceMeshComplex_2InTriangulation_3`, 
a data structure able to represent a two dimensional 
complex embedded in a three dimensional triangulation. 
The argument `c2t3` of type `SurfaceMeshC2T3`, passed by reference 
to the surface mesh generator, 
is used to maintain the current approximating mesh and it stores 
the final mesh at the end of the procedure. 
The type `SurfaceMeshC2T3` is in particular required to 
provide a type `SurfaceMeshC2T3::Triangulation_3` 
for the three dimensional triangulation 
embedding the surface mesh. 
The vertex and cell base classes of the triangulation 
`SurfaceMeshC2T3::Triangulation_3` are required 
to be models of the concepts 
`SurfaceMeshVertexBase_3` and 
`SurfaceMeshCellBase_3` respectively. 

\tparam Surface stands for the surface type. 
This type has to be a model of the concept `Surface_3`. 
<br>
The knowledge on the surface, required by the surface mesh generator 
is encapsulated in a 
traits class. Actually, the mesh generator accesses the surface to be meshed 
through this traits class only. 
The traits class is required to be a model 
of the concept `SurfaceMeshTraits_3`. 

\tparam FacetsCriteria has to be a model 
of the concept `SurfaceMeshFacetsCriteria_3`. 
The argument of type `FacetsCriteria` passed to the surface 
mesh generator specifies the size and shape requirements 
on the output surface mesh. 

\tparam Tag 
is a tag whose type affects the behavior of the 
meshing algorithm. The function `make_surface_mesh()` has specialized versions 
for the following tag types: 
<br>
- `Manifold_tag`: the output mesh is guaranteed to be a manifold 
surface without boundary. 
- `Manifold_with_boundary_tag`: the output mesh is guaranteed to be 
manifold but may have boundaries. 
- `Non_manifold_tag`: the algorithm relies on the given criteria and 
guarantees nothing else. 

The Delaunay refinement 
process is started with an initial set of points which is the union 
of two sets: the 
set of vertices in the initial triangulation pointed to by the 
`c2t3` argument and a set of 
points provided by the traits class. 
The optional parameter `initial_number_of_points` 
allows to monitor the number of points in this second set. 
(This parameter is passed to the `operator()` of 
the constructor object `Construct_initial_points` 
in the traits class.) 
The meshing algorithm requires that the initial set of points 
includes at least one point 
on each connected components of the surface to be meshed. 


\sa `SurfaceMeshComplex_2InTriangulation_3` 
\sa `SurfaceMeshCellBase_3` 
\sa `SurfaceMeshVertexBase_3` 
\sa `Surface_3` 
\sa `SurfaceMeshFacetsCriteria_3` 
\sa `Surface_mesh_default_triangulation_3` 

*/

/// @{

/*!
In the first overloaded version of 
of `make_surface_mesh()`, the surface type is given 
as template parameter (`Surface`) and the `surface` 
to be meshed is passed as parameter to the mesh generator. 
In that case the surface mesh generator traits type 
is automatically generated form the surface type 
by an auxiliary class called the `Surface_mesh_traits_generator_3`.


The first overloaded version can be used 
whenever the surface type either provides a nested type 
`Surface::Surface_mesher_traits_3` 
that is a model of `SurfaceMeshTraits_3` 
or is a surface type for which a specialization 
of the traits generator `Surface_mesh_traits_generator_3<Surface>` 
is provided. 
Currently, the library provides partial specializations 
of `Surface_mesh_traits_generator_3<Surface>` 
for implicit surfaces (`Implicit_surface_3<Traits, Function>`) and 
gray level images (`Gray_level_image_3<FT, Point>`). 

 */  
template <class SurfaceMeshC2T3,
class Surface,
class FacetsCriteria,
class Tag >
void make_surface_mesh(SurfaceMeshC2T3& c2t3,
Surface surface,
FacetsCriteria criteria,
Tag tag,
int initial_number_of_points = 20) ;

/*!
In the second overloaded version of `make_surface_mesh()`, 
the surface mesh generator traits type is provided 
by the template parameter `SurfaceMeshTraits_3` 
and the surface type is obtained from this traits type. 
Both the surface and the traits 
are passed to the mesh generator as arguments. 
 */
template <class SurfaceMeshC2T3,
class SurfaceMeshTraits,
class FacetsCriteria,
class Tag >
void make_surface_mesh(SurfaceMeshC2T3& c2t3,
SurfaceMeshTraits::Surface_3 surface,
SurfaceMeshTraits traits,
FacetsCriteria criteria,
Tag tag,
int initial_number_of_points = 20 );

/// @}
} /* namespace CGAL */


namespace CGAL {

/*!
\ingroup PkgSurfaceMesher3TagClasses

The class `Manifold_tag` is a tag class used to monitor the 
surface meshing algorithm. When instantiated with the tag 
`Manifold_tag` the function template 
`make_surface_mesh()` 
ensures that the output mesh is a manifold surface 
without boundary. 

\sa `make_surface_mesh()` 
\sa `Manifold_with_boundary_tag` 
\sa `Non_manifold_tag` 

*/

class Manifold_tag {
public:

/// @}

}; /* end Manifold_tag */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSurfaceMesher3TagClasses

The class `Manifold_with_boundary_tag` is a tag class used to monitor the 
surface meshing algorithm. When instantiated with the tag 
`Manifold_with_boundary_tag`, the function template 
`make_surface_mesh()` 
ensures that the output mesh is a manifold surface 
but it may have boundaries. 

\sa `make_surface_mesh()` 
\sa `Manifold_tag` 
\sa `Non_manifold_tag` 

*/

class Manifold_with_boundary_tag {
public:

/// @}

}; /* end Manifold_with_boundary_tag */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSurfaceMesher3TagClasses

The class `Non_manifold_tag` is a tag class used to monitor the 
surface meshing algorithm. When instantiated with the tag 
`Non_manifold_tag` the function template 
`make_surface_mesh()` 
does not ensure that the output mesh is a manifold surface. 
The manifold property of output mesh 
may nevertheless result from the choice of 
appropriate meshing criteria. 

\sa `make_surface_mesh()` 
\sa `Manifold_tag` 
\sa `Manifold_with_boundary_tag` 

*/

class Non_manifold_tag {
public:

/// @}

}; /* end Non_manifold_tag */
} /* end namespace CGAL */
