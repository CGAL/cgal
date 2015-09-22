namespace CGAL {

/*!
\ingroup kernel_classes2

An object `r` of the data type `Iso_rectangle_2` is a 
rectangle in the Euclidean plane \f$ \E^2\f$ with sides parallel to the \f$ x\f$ and 
\f$ y\f$ axis of the coordinate system. 

Although they are represented in a canonical form by only two 
vertices, namely the lower left and the upper right vertex, we provide 
functions for "accessing" the other vertices as well. The vertices 
are returned in counterclockwise order. 

Iso-oriented rectangles and bounding boxes are quite similar. The 
difference however is that bounding boxes have always double coordinates, 
whereas the coordinate type of an iso-oriented rectangle is chosen by 
the user. 

\sa `Kernel::IsoRectangle_2` 

*/
template< typename Kernel >
class Iso_rectangle_2 {
public:

/// \name Creation 
/// @{

/*!
introduces an iso-oriented rectangle `r` with diagonal 
opposite vertices `p` and `q`. Note that the object is 
brought in the canonical form. 
*/ 
Iso_rectangle_2(const Point_2<Kernel> &p, 
const Point_2<Kernel> &q); 

/*!
introduces an iso-oriented rectangle `r` with diagonal 
opposite vertices `p` and `q`. The `int` argument value 
is only used to distinguish the two overloaded functions. 
\pre `p.x()<=q.x()` and `p.y()<=q.y()`. 
*/ 
Iso_rectangle_2(const Point_2<Kernel> &p, 
const Point_2<Kernel> &q, 
int); 

/*!
introduces an iso-oriented rectangle `r` whose 
minimal \f$ x\f$ coordinate is the one of `left`, the 
maximal \f$ x\f$ coordinate is the one of `right`, the 
minimal \f$ y\f$ coordinate is the one of `bottom`, the 
maximal \f$ y\f$ coordinate is the one of `top`. 
*/ 
Iso_rectangle_2(const Point_2<Kernel> &left, 
const Point_2<Kernel> &right, 
const Point_2<Kernel> &bottom, 
const Point_2<Kernel> &top); 

/*!
introduces an iso-oriented rectangle `r` with diagonal 
opposite vertices (`min_hx/hw`, `min_hy/hw`) and 
(`max_hx/hw`, `max_hy/hw`). 
\pre `hw` \f$ \neq\f$ 0. 
*/ 
Iso_rectangle_2(const Kernel::RT& min_hx, const Kernel::RT& min_hy, 
const Kernel::RT& max_hx, const Kernel::RT& max_hy, 
const Kernel::RT& hw = RT(1)); 

/*!
If `Kernel::RT` is constructible from double,
introduces an iso-oriented rectangle from `bbox`.
*/
Iso_rectangle_2(const Bbox_2& bbox);

/// @} 

/// \name Operations 
/// @{

/*!
Test for equality: two iso-oriented rectangles are equal, iff their 
lower left and their upper right vertices are equal. 
*/ 
bool operator==(const Iso_rectangle_2<Kernel> &r2) const; 

/*!
Test for inequality. 
*/ 
bool operator!=(const Iso_rectangle_2<Kernel> &r2) const; 

/*!
returns the i'th vertex modulo 4 of `r` in counterclockwise order, 
starting with the lower left vertex. 
*/ 
Point_2<Kernel> vertex(int i) const; 

/*!
returns `vertex(i)`. 
*/ 
Point_2<Kernel> operator[](int i) const; 

/*!
returns the lower left vertex of `r` (= `vertex(0)`). 
*/ 
Point_2<Kernel> min() const; 

/*!
returns the upper right vertex of `r` (= `vertex(2)`). 
*/ 
Point_2<Kernel> max() const; 

/*!
returns the \f$ x\f$ coordinate of lower left vertex of `r`. 
*/ 
Kernel::FT xmin() const; 

/*!
returns the \f$ y\f$ coordinate of lower left vertex of `r`. 
*/ 
Kernel::FT ymin() const; 

/*!
returns the \f$ x\f$ coordinate of upper right vertex of `r`. 
*/ 
Kernel::FT xmax() const; 

/*!
returns the \f$ y\f$ coordinate of upper right vertex of `r`. 
*/ 
Kernel::FT ymax() const; 

/*!
returns the `i`'th %Cartesian coordinate of the 
lower left vertex of `r`. 
\pre \f$ 0 \leq i \leq1\f$. 
*/ 
Kernel::FT min_coord(int i) const; 

/*!
returns the `i`'th %Cartesian coordinate of the 
upper right vertex of `r`. 
\pre \f$ 0 \leq i \leq1\f$. 
*/ 
Kernel::FT max_coord(int i) const; 

/// @} 

/// \name Predicates 
/// @{

/*!
`r` is degenerate, if all vertices 
are collinear. 
*/ 
bool is_degenerate() const; 

/*!
returns either \ref ON_UNBOUNDED_SIDE, 
\ref ON_BOUNDED_SIDE, or the constant 
\ref ON_BOUNDARY, 
depending on where point `p` is. 
*/ 
Bounded_side bounded_side(const Point_2<Kernel> &p) const; 

/*!

*/ 
bool has_on_boundary(const Point_2<Kernel> &p) const; 

/*!

*/ 
bool has_on_bounded_side(const Point_2<Kernel> &p) const; 

/*!

*/ 
bool has_on_unbounded_side(const Point_2<Kernel> &p) const; 

/// @} 

/// \name Miscellaneous 
/// @{

/*!
returns the area of `r`. 
*/ 
Kernel::FT area() const; 

/*!
returns a bounding box containing `r`. 
*/ 
Bbox_2 bbox() const;

/*!
returns the iso-oriented rectangle obtained by applying `t` on 
the lower left and the upper right corner of `r`. 
\pre The angle at a rotation must be a multiple of \f$ \pi/2\f$, otherwise the resulting rectangle does not have the same side length. Note that rotating about an arbitrary angle can even result in a degenerate iso-oriented rectangle. 
*/ 
Iso_rectangle_2<Kernel> transform(const Aff_transformation_2<Kernel> &t) const; 

/// @}

}; /* end Iso_rectangle_2 */
} /* end namespace CGAL */
