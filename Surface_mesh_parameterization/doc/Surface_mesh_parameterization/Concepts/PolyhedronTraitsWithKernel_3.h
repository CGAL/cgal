
/*!
\ingroup PkgSurfaceParameterizationConcepts
\cgalConcept

\anchor pagePolyTraitsRef 

Required types for the `PolyhedronTraitsWithKernel_3` concept. This 
geometric traits concept is used in the polyhedral surface data 
structure `CGAL::Polyhedron_3<Traits>`. This concept provides 
additional requirements to the `PolyhedronTraits_3` concept 
required by `CGAL::Polyhedron_3<Traits>` used within the class 
`CGAL::Parameterization_polyhedron_adaptor_3<Polyhedron_3_>`. 

\cgalRefines `PolyhedronTraits_3` 

\cgalHasModel `CGAL::Polyhedron_traits_3<Kernel>`
\cgalHasModel `CGAL::Polyhedron_traits_with_normals_3<Kernel>`
\cgalHasModel All models of the `CGAL::Kernel` concept such as `Simple_cartesian<FieldNumberType>`. 

\sa `CGAL::Parameterization_polyhedron_adaptor_3<Polyhedron_3_>`
\sa `CGAL::Polyhedron_3<Traits>`

*/

class PolyhedronTraitsWithKernel_3 {
public:

/// \name Types 
/// @{

/*!
a kernel type providing: a field type `FT`; 2D and 3D point types `Point_2` 
and `Point_3`; 2D and 3D vector types `Vector_2` and `Vector_3`. 
*/ 
typedef unspecified_type Kernel; 

/// @}

}; /* end PolyhedronTraitsWithKernel_3 */

