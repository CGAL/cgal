/*!
\ingroup PkgSpatialSearchingDConcepts
\cgalConcept

The concept `RangeSearchTraits` defines the requirements for the template 
parameter of the search classes. This concept also defines requirements to 
range search queries in a model of `SpatialTree`. 

\cgalRefines `SearchTraits` 

\cgalHasModel `CGAL::Cartesian_d<FT>` 
\cgalHasModel `CGAL::Homogeneous_d<RT>` 
\cgalHasModel `CGAL::Search_traits_2<Kernel>` 
\cgalHasModel `CGAL::Search_traits_3<Kernel>` 

\sa `SearchTraits` 
\sa `CGAL::Search_traits_adapter<Key,PointPropertyMap,BaseTraits>` 

*/

class RangeSearchTraits {
public:

/// \name Types 
/// @{

/*!
Iso box type, which is only needed for range search queries. 
*/ 
typedef unspecified_type Iso_box_d; 

/*!
Sphere type, which is only needed for range search queries. 
*/ 
typedef unspecified_type Sphere_d; 

/*!
Functor with operator to construct the iso box from two points. 
*/ 
typedef unspecified_type Construct_iso_box_d; 

/*!
Functor with operator to construct 
the center of an object of type `Sphere_d`. 
*/ 
typedef unspecified_type Construct_center_d; 

/*!
Functor with operator to compute 
the squared radius of a an object of type `Sphere_d`. 
*/ 
typedef unspecified_type Compute_squared_radius_d;

/*!
Functor with operator to construct 
the vertex with lexicographically smallest coordinates of an object of type `Iso_box_d`. 
*/ 
typedef unspecified_type Construct_min_vertex_d; 

/*!
Functor with operator to construct 
the vertex with lexicographically largest coordinates of an object of type `Iso_box_d`. 
*/ 
typedef unspecified_type Construct_max_vertex_d; 

/// @}

}; /* end RangeSearchTraits */
