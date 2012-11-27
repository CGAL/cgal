/*!
\ingroup PkgGeneratorsConcepts
\cgalConcept

The concept `RandomConvexSetTraits_2` describes the requirements of the traits 
class for the function `random_convex_set_2()`. 

\cgalHasModel `CGAL::Random_convex_set_traits_2<Kernel>`

*/

class RandomConvexSetTraits_2 {
public:

/// \name Types 
/// @{

/*! 
point class. 
*/ 
typedef Hidden_type Point_2; 

/*! 
class used for doing computations on point and 
vector coordinates (has to fulfill field type requirements). 
*/ 
typedef Hidden_type FT; 

/*! 
AdaptableBinaryFunction class: 
`Point_2` \f$ \times\f$ `Point_2` \f$ \rightarrow\f$ 
`Point_2`. It returns the point that results from adding 
the vectors corresponding to both arguments. 
*/ 
typedef Hidden_type Sum; 

/*! 
AdaptableBinaryFunction class: 
`Point_2` \f$ \times\f$ `FT` \f$ \rightarrow\f$ 
`Point_2`. `Scale(p,k)` returns the point that 
results from scaling the vector corresponding to `p` by a 
factor of `k`. 
*/ 
typedef Hidden_type Scale; 

/*! 
AdaptableUnaryFunction class: 
`Point_2` \f$ \rightarrow\f$ `FT`. `Max_coordinate(p)` 
returns the coordinate of `p` with largest absolute value. 
*/ 
typedef Hidden_type Max_coordinate; 

/*! 
AdaptableBinaryFunction class: 
`Point_2` \f$ \times\f$ `Point_2` \f$ \rightarrow\f$ 
`bool`. It returns `true`, iff the angle of the 
direction corresponding to the first argument with respect to 
the positive \f$ x\f$-axis is less than the angle of the direction 
corresponding to the second argument. 
*/ 
typedef Hidden_type Angle_less; 

/// @} 

/// \name Operations 
/// @{

/*! 
return origin (neutral 
element for the `Sum` operation). 
*/ 
Point_2 origin() const; 

/// @}

}; /* end RandomConvexSetTraits_2 */
