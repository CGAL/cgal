
namespace CGAL {

/*!
\ingroup PkgKernel23

An object `b` of the class `Bbox_3` is a bounding 
box in the three-dimensional Euclidean space \f$ \E^3\f$. 

\sa `CGAL::Bbox_2` 
\sa `CGAL::do_overlap` 

*/

class Bbox_3 {
public:

/// \name Creation 
/// @{

/*! 
introduces a bounding box `b` with lexicographically 
smallest corner point at `(xmin, ymin, zmin)` 
lexicographically largest corner point at 
`(xmax, ymax, zmax)`. 
*/ 
Bbox_3(double x_min, double y_min, double z_min, 
double x_max, double y_max, double z_max); 

/// @} 

/// \name Operations 
/// @{

/*! 
Test for equality. 
*/ 
bool operator==(const Bbox_3 &c) const; 

/*! 
Test for inequality. 
*/ 
bool operator!=(const Bbox_3 &q) const; 

/*! 
Returns 3. 
*/ 
int dimension() const; 

/*! 

*/ 
double xmin() const; 

/*! 

*/ 
double ymin() const; 

/*! 

*/ 
double zmin() const; 

/*! 

*/ 
double xmax() const; 

/*! 

*/ 
double ymax() const; 

/*! 

*/ 
double zmax() const; 

/*! 
Returns `xmin()` if `i==0` or `ymin()` if `i==1` 
or `zmin()` if `i==2`. 
\pre i<=0 and i<=2 
*/ 
double min(int i) const; 

/*! 
Returns `xmax()` if `i==0` or `ymax()` if `i==1` 
or `zmax()` if `i==2`. 
\pre i==0 and i<=2 
*/ 
double max(int i) const; 

/*! 
returns a bounding box of `b` and `c`. 
*/ 
Bbox_3 operator+(const Bbox_3 &c) const; 

/// @}

}; /* end Bbox_3 */

/*!
returns `true` iff `bb1` and `bb2` overlap, i.e., iff their
intersection is non-empty.

\relates Bbox_3
*/
bool do_overlap(const Bbox_3 &bb1, const Bbox_3 &bb2);

} /* end namespace CGAL */
