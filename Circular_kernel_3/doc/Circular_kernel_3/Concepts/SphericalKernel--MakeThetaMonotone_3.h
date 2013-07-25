
/*!
\ingroup PkgSphericalKernel3GeometricConcepts
\cgalConcept


\sa `SphericalKernel::IsThetaMonotone_3`

*/

class SphericalKernel::MakeThetaMonotone_3 {
public:

/// \name Operations
/// A model of this concept must provide: 
/// @{

/*!

Copies in the output iterator the results of the split of arc `a` at the \f$ \theta\f$-extremal 
point(s) of its supporting circle relatively to the context sphere used by the function `SphericalKernel::make_theta_monotone_3_object` 
(Refer to Section \ref sectionSKobjects for the definition of these points.) 
The output iterator may contain no arc (if the supporting circle is a bipolar circle), 
one arc (if `a` is already \f$ \theta\f$-monotone), two arcs (if only one \f$ \theta\f$-extremal point is on `a`), or 
three arcs (if two \f$ \theta\f$-extremal points are on `a`). 
\pre `a` lies on the context sphere used by the function `SphericalKernel::make_theta_monotone_3_object`, and the supporting circle of `a` is not bipolar. 

*/ 
template<class OutputIterator> 
OutputIterator operator() 
(const SphericalKernel::Circular_arc_3 &a,OutputIterator res); 

/*!
Copies in the output iterator the results of the split of circle `c` at its \f$ \theta\f$-extremal 
point(s) relatively to the context sphere used by the function `SphericalKernel::make_theta_monotone_3_object`. 
(Refer to Section \ref sectionSKobjects for the definition of these points.) 
The output iterator may contain no arc (if the circle is bipolar), 
one arc (if the circle is polar or threaded), or two arcs (if the circle is normal). 

The source and target are such that 
the circular arc is the set of points of the circle that lie between the source 
and the target when traversing the circle counterclockwise 
seen from the positive side of the plane of the circle. 

In this definition, we say that a normal vector \f$ (a,b,c)\f$ is <I>positive</I> if 
\f$ (a>0) || (a==0) \&\& (b>0) || (a==0)\&\&(b==0)\&\&(c>0)\f$. 

For a threaded circle, the arc returned the one built using the full circle. 

For a polar circle, the arc returned is the full circle, the source and target correspond to the pole the circle goes through. 
\pre `c` lies on the context sphere used by the function `SphericalKernel::make_theta_monotone_3_object`, and `c` is not bipolar. 

*/ 
template<class OutputIterator> 
OutputIterator operator() 
(const SphericalKernel::Circle_3 &c,OutputIterator res); 

/// @}

}; /* end SphericalKernel::MakeThetaMonotone_3 */

