/*!
\ingroup PkgSpatialSearchingDConcepts
\cgalconcept

The concept `RangeSearchTraits` defines the requirements for the template 
parameter of the search classes. This concept also defines requirements to 
range search queries in a model of `SpatialTree`. 

\refines ::SearchTraits 

\hasModel CGAL::Cartesian_d<FT> 
\hasModel CGAL::Homogeneous_d<RT> 
\hasModel CGAL::Search_traits_2<Kernel> 
\hasModel CGAL::Search_traits_3<Kernel> 

\sa `SearchTraits` 
\sa `CGAL::Search_traits_adapter<Key,PointPropertyMap,BaseTraits>` 

*/

class RangeSearchTraits {
public:

/// \name Types 
/// @{

/*! 
Iso box type. It is only needed for range search queries. 
*/ 
typedef Hidden_type Iso_box_d; 

/*! 
Sphere type. It is only needed for range search queries. 
*/ 
typedef Hidden_type Sphere_d; 

/*! 
Functor with operator to construct the iso box from two points. 
*/ 
typedef Hidden_type Construct_iso_box_d; 

/*! 
Functor with operator to construct 
the center of an object of type `Sphere_d`. 
*/ 
typedef Hidden_type Construct_center_d; 

/*! 
Functor with operator to compute 
the squared radius of a an object of type `Sphere_d`. 
*/ 
typedef Hidden_type Construct_squared_radius_d; 

/*! 
Functor with operator to construct 
the vertex with lexicographically smallest coordinates of an object of type `Iso_box_d`. 
*/ 
typedef Hidden_type Construct_min_vertex_d; 

/*! 
Functor with operator to construct 
the vertex with lexicographically largest coordinates of an object of type `Iso_box_d`. 
*/ 
typedef Hidden_type Construct_max_vertex_d; 

/// @}

}; /* end RangeSearchTraits */
