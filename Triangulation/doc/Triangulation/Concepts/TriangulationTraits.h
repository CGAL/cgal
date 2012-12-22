
/*!
\ingroup PkgTriangulationsConcepts
\cgalconcept

The concept `TriangulationTraits` describes the various types and functions that a class 
must provide as the first parameter (`TriangulationTraits`) to the class template 
`Triangulation<TriangulationTraits, TriangulationDataStructure>`. It brings the geometric ingredient to the 
definition of a triangulation, while the combinatorial ingredient is brought by 
the second template parameter, `TriangulationDataStructure`. 

Inserting a range of points in a triangulation is optimized using 
spatial sorting, thus besides the requirements below, 
a class provided as `TriangulationTraits` should also satisfy the concept 
`SpatialSortingTraits_d`. 

\refines ::SpatialSortingTraits_d 
CONVERRORRefines: If a range of points is inserted, the 
traits must refine \refines ::SpatialSortingTraits_d, This is not needed 
CONVERRORRefines: if the points are inserted one by one. 

CONVERRORHasModels: `CGAL::Cartesian_d<FT, Dim, LA>`, 
CONVERRORHasModels: `CGAL::????<K>` (recommended). 

\sa `DelaunayTriangulationTraits` 
\sa `Triangulation` 

*/

class TriangulationTraits {
public:

/// \name Types 
CONVERROR Check if this needs to be spread\n/// In the \f$ D\f$-dimensional oriented space, a \f$ k-1\f$ dimensional subspace (flat) define by \f$ k\f$ points can be oriented in two different ways. Choosing the orientation of any simplex defined by \f$ k\f$ points fix the orientation of all other simplices. To be able to orient lower dimensional flats, we use the following classes:
/// @{

/*! 
A type representing the dimension of the underlying space. it can be static 
(`Maximal_dimension`=`CGAL::``Dimension_tag<int dim>`) or 
dynamic (`Maximal_dimension`=`CGAL::``Dynamic_dimension_tag`). 
This dimension must match the dimension of the predicate 
`Orientation_d` but not necessarily the one of `Point_d`. 

*/ 
typedef Hidden_type Dimension; 

/*! 
A type representing a point in Euclidean space. 
*/ 
typedef Hidden_type Point_d; 

/*! 
Functor returning the dimension of a `Point_d`. 
Must provide 
`int operator()(Point_d p)` returning the dimension of \f$ p\f$. 

*/ 
typedef Hidden_type Point_dimension_d; 

/*! 
A predicate object that must provide the 
templated operator 
`template<typename ForwardIterator> Orientation operator()(ForwardIterator start, ForwardIterator end)`. 
The operator returns 
`CGAL::POSITIVE`, `CGAL::NEGATIVE` or `CGAL::COPLANAR` depending on 
the orientation of the simplex defined by the points in the range `[start, end)`. 
\pre `std::distance(start,end)=D+1`, where 
`Point_dimension_d(*it)` is \f$ D\f$ for all `it` in `[start,end)`. 

*/ 
typedef Hidden_type Orientation_d; 

/*! 
A predicate object that must provide 
the templated operator 
`template<typename ForwardIterator> bool operator()(ForwardIterator start, ForwardIterator end, const Point_d & p)`. 
The operator returns `true` if and only if point `p` is 
contained in the affine space spanned by the points in the range `[start, end)`. That affine space is also called the <I>affine hull</I> of the points 
in the range. 
\pre The \f$ k\f$ points in the range 
must be affinely independent. 
`Point_dimension_d(*it)` is \f$ D\f$ for all `it` in 
`[start,end)`, for some \f$ D\f$. 
\f$ 2\leq k\leq D\f$. 

*/ 
typedef Hidden_type Contained_in_affine_hull_d; 

/*! 

A type representing an orientation of an affine subspace of 
dimension \f$ k\f$ strictly smaller than the maximal dimension. 

*/ 
typedef Hidden_type Flat_orientation_d; 

/*! 

A construction object that must 
provide the templated operator 
`template<typename ForwardIterator> Flat_orientation_d operator()(ForwardIterator start, ForwardIterator end)`. 

The flat spanned by the points in 
the range `R=[start, end)` can be oriented in two different ways, 
the operator 
returns an object that allow to orient that flat so that `R=[start, end)` 
defines a positive simplex. 
\pre The \f$ k\f$ points in the range 
must be affinely independent. 
`Point_dimension_d(*it)` is \f$ D\f$ for all `it` in `R` for 
some \f$ D\f$. 
\f$ 2\leq k\leq D\f$. 

*/ 
typedef Hidden_type Construct_flat_orientation_d; 

/*! 

A predicate object that must provide the 
templated operator 
`template<typename ForwardIterator> Orientation operator()(Flat_orientation_d orient,ForwardIterator start, ForwardIterator end)`. 

The operator returns 
`CGAL::POSITIVE`, `CGAL::NEGATIVE` or `CGAL::COPLANAR` depending on 
the orientation of the simplex defined by the points in the range `[start, end)`. 
The points are supposed to belong to the lower dimensional flat 
whose orientation is given by `orient`. 
\pre `std::distance(start,end)=k` where \f$ k\f$ is the number of 
points 
used to construct `orient`. 
`Point_dimension_d(*it)` is \f$ D\f$ for all `it` in 
`[start,end)` where \f$ D\f$ is the dimension of the points used to 
construct `orient`. 
\f$ 2\leq k\leq D\f$. 

*/ 
typedef Hidden_type In_flat_orientation_d; 

/*! 
A predicate object that must 
provide the operator 
`Comparison_result operator()(const Point_d & p, const Point_d & q)`. 
The operator returns `SMALLER` if `p` is 
lexicographically smaller than point `q`, `EQUAL` if both points are 
the same and `LARGER` otherwise. 
*/ 
typedef Hidden_type Compare_lexicographically_d; 

/// @} 

/// \name Creation 
/// @{

/*! 
The default constructor. 
*/ 
TriangulationTraits(); 

/// @} 

/// \name Operations 
CONVERROR Check if this needs to be spread\n/// The following methods permit access to the traits class's predicates:
/// @{

/*! 

*/ 
Orientation_d orientation_d_object() const; 

/*! 

*/ 
Contained_in_affine_hull_d contained_in_affine_hull_d_object() 
const; 

/*! 

*/ 
Construct_flat_orientation_d construct_flat_orientation_d_object() const; 

/*! 

*/ 
In_flat_orientation_d in_flat_orientation_d_object() const; 

/*! 

*/ 
Compare_lexicographically_d compare_lexicographically_d_object() 
const; 

/// @}

}; /* end TriangulationTraits */

