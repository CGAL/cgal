
/*!
\ingroup PkgTriangulationsConcepts
\cgalConcept

This concept describes the geometric types and predicates required to build
a triangulation. It corresponds to the first template parameter of the class 
`CGAL::Triangulation<TriangulationTraits, TriangulationDataStructure>`.

\cgalRefines `SpatialSortingTraits_d`

If a range of points is inserted, the 
traits must refine `SpatialSortingTraits_d`. The insertion is 
then optimized using spatial sorting. 
This is not required if the points are inserted one by one.

\cgalHasModel `CGAL::Epick_d<Dim>`

\sa `DelaunayTriangulationTraits`
*/

class TriangulationTraits {
public:

/// \name Types 
/// @{

/*!
A type representing the dimension of the predicates 
(but not necessarily the one of `Point_d`). If \f$ n \f$ is the number of
points required by the `Orientation_d` predicate, then 
`Dimension` \f$ = n - 1\f$.
It can be static (`Dimension`=`CGAL::``Dimension_tag<int dim>`) or 
dynamic (`Dimension`=`CGAL::``Dynamic_dimension_tag`).
*/ 
typedef unspecified_type Dimension;

/*! 
A type representing a point in Euclidean space. It must be 
`DefaultConstructible`, `CopyConstructible` and `Assignable`.
*/ 
typedef unspecified_type Point_d; 

/*! 
A predicate object that must provide the 
templated operator 
`template<typename ForwardIterator> Orientation operator()(ForwardIterator start, ForwardIterator end)`.
The operator returns the orientation of the simplex defined by the points 
in the range `[start, end)`; the value can be 
`CGAL::POSITIVE`, `CGAL::NEGATIVE` or `CGAL::COPLANAR`.
\pre If `Dimension`=`CGAL::``Dimension_tag<D>`, then `std::distance(start,end)=D+1`.
*/ 
typedef unspecified_type Orientation_d; 

/*! 
A predicate object that must provide 
the templated operator 
`template<typename ForwardIterator> bool operator()(ForwardIterator start, ForwardIterator end, const Point_d & p)`. 
The operator returns `true` if and only if point `p` is 
contained in the affine space spanned by the points in the range `[start, end)`. That affine space is also called the <I>affine hull</I> of the points 
in the range.
\pre If `Dimension`=`CGAL::``Dimension_tag<D>`, 
then `std::distance(start,end)=D+1`.
The points in the range 
must be affinely independent. Note that in the CGAL kernels, this predicate
works also with affinely dependent points.
\f$ 2\leq k\leq D\f$. 

*/ 
typedef unspecified_type Contained_in_affine_hull_d; 

/// @}

/// \name
/// In the \f$ D\f$-dimensional oriented space, a \f$ k-1\f$
/// dimensional subspace (flat) defined by \f$ k\f$ points can be
/// oriented in two different ways. Choosing the orientation of any
/// simplex defined by \f$ k\f$ points in a flat fixes the orientation of
/// the flat and therefore the orientation of all other simplices in this flat.
/// To be able to orient lower dimensional flats, we
/// use the following classes:
/// @{

/*!
A type representing an orientation of an affine subspace of 
dimension \f$ k\f$ strictly smaller than the dimension of the traits. 
*/ 
typedef unspecified_type Flat_orientation_d; 

/*!
A construction object that must 
provide the templated operator 
`template<typename ForwardIterator> Flat_orientation_d operator()(ForwardIterator start, ForwardIterator end)`. 

The flat spanned by the points in 
the range `R=[start, end)` can be oriented in two different ways, 
the operator 
returns an object that allow to orient that flat so that `R=[start, end)` 
defines a positive simplex.
\pre If `Dimension`=`CGAL::``Dimension_tag<D>`, 
then `std::distance(start,end)=D+1`.
The points in range
`[start,end)` must be affinely independent.
\f$ 2\leq k\leq D\f$. 
*/ 
typedef unspecified_type Construct_flat_orientation_d; 

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
points used to construct `orient`. 
\f$ 2\leq k\leq D\f$. 
*/ 
typedef unspecified_type In_flat_orientation_d; 

/*!
A predicate object that must 
provide the operator 
`Comparison_result operator()(const Point_d & p, const Point_d & q)`. 
The operator returns `SMALLER` if `p` is 
lexicographically smaller than point `q`, `EQUAL` if both points are 
the same and `LARGER` otherwise. 
*/ 
typedef unspecified_type Compare_lexicographically_d; 

/// @} 

/// \name Creation 
/// @{

/*! 
The default constructor. 
*/ 
TriangulationTraits(); 

/// @} 

/// \name Operations 
/// The following methods permit access to the traits class's predicates:
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

