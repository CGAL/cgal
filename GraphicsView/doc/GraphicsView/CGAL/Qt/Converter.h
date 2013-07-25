namespace CGAL {
namespace Qt {
/*!
\ingroup PkgGraphicsViewMiscClasses

An object of type `Converter` converts 2D \cgal `Kernel` objects to their equivalent 
objects in <I>Qt</I>, and the other way round. Note that some objects have no equivalent. 
For example the `CGAL::Circle_2<K>` cannot be converted to something in Qt, 
and the unbounded objects `CGAL::Ray_2<K>` and `CGAL::Line_2<K>` are clipped. 
Note also that \cgal and <I>Qt</I> sometimes also use the same word for different things. 
For example <I>line</I> denotes an unbounded line in \cgal, wheras it denotes a bounded 
segment in <I>Qt</I>. 

\tparam K must be a model of `Kernel`.

*/
template< typename K >
class Converter {
public:

/// \name Creation 
/// @{

/*!
The clipping rect is used for unbounded \cgal objects. 
*/ 
Converter<K>(QRectF clippingRect); 

/// @} 

/// \name From CGAL to Qt 
/// @{

/*!
Converts a bounding box. 
*/ 
QRectF operator()(Bbox_2); 

/*!
Converts a point. 
*/ 
QPointF operator()(K::Point_2); 

/*!
Converts a segment. 
*/ 
QLineF operator()(K::Segment_2); 

/*!
Clips the ray and converts the resulting segment. 
*/ 
QLineF operator()(K::Ray_2); 

/*!
Clips the line and converts the resulting segment. 
*/ 
QLineF operator()(K::Line_2); 

/*!
Converts a triangle. 
*/ 
QPolygonF operator()(K::Triangle_2); 

/*!
Converts an iso rectangle. 
*/ 
QRectF operator()(K::Iso_rectangle_2); 

/*!
Converts a list of points to a polygon. 
*/ 
QPolygonF operator()(std::list<K::Point_2>); 

/// @} 

/// \name From Qt to CGAL 
/// @{

/*!
Converts a point. 
*/ 
K::Point_2 operator()(QPointF); 

/*!
Converts a segment. 
*/ 
K::Segment_2 operator()(QLineF); 

/*!
Converts an iso rectangle. 
*/ 
K::Iso_rectangle_2 operator()(QRectF); 

/*!
Converts a polygon to a list of points. 
*/ 
std::list<K::Point_2> operator()(QPolygonF); 

/// @}

}; /* end Converter */
} /* end namespace Qt */
} /* end namespace CGAL */
