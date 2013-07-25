namespace CGAL {
namespace Qt {
/*!
\ingroup PkgGraphicsViewMiscClasses

An object of type `PainterOstream` provides output operators for `Kernel` objects. 
As \cgal has unbounded objects the `PainterOstream` must be constructed with a clipping rectangle. 
This is typically the visible area of the widget where the unbounded object is displayed. 

*/
template< typename K >
class PainterOstream {
public:

/// \name Creation 
/// @{

/*!
The clipping rect is used for unbounded `Kernel` objects. 
*/ 
PainterOstream<K>(QPainter* qp, QRectF clippingRect); 

/// @} 

/// \name Operations 
/// @{

/*!
Draws a point. 
*/ 
PainterOstream<K> operator<<(K::Point_2); 

/*!
Draws a segment. 
*/ 
PainterOstream<K> operator<<(K::Segment_2); 

/*!
Draws a clipped ray. 
*/ 
PainterOstream<K> operator<<(K::Ray_2); 

/*!
Draws a clipped line. 
*/ 
PainterOstream<K> operator<<(K::Line_2); 

/*!
Draws a triangle. 
*/ 
PainterOstream<K> operator<<(K::Triangle_2); 

/*!
Draws an iso rectangle. 
*/ 
PainterOstream<K> operator<<(K::Iso_rectangle_2); 

/*!
Draws a circle. 
*/ 
PainterOstream<K> operator<<(K::Circle_2); 

/*!
Draws a circular arc. 
*/ 
PainterOstream<K> operator<<(K::Circular_arc_2); 

/*!
Draws a polyline. In order to close it the first and last point must be equal. 
*/ 
PainterOstream<K> operator<<(std::list<K::Point_2>); 

/*!
Draws an iso rectangle. 
*/ 
PainterOstream<K> operator<<(Bbox_2); 

/*!
Sets the pen used in the next paint operations. 
*/ 
PainterOstream<K> operator<<(QPen); 

/*!
Sets the brush used in the next paint operations. 
*/ 
PainterOstream<K> operator<<(QBrush); 

/// @}

}; /* end PainterOstream */
} /* end namespace Qt */
} /* end namespace CGAL */
