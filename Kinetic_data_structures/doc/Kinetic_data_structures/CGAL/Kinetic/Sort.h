
namespace CGAL {
namespace Kinetic {

/*!
\ingroup PkgKdsSorting

The class `Kinetic::Sort` maintains a sorted list of objects. It is the 
simplest kinetic data structure provided and is a good place to start 
when looking at the basics of implementing a kinetic data 
structure. 

The `Kinetic::SortVisitor` can be used to monitor what is happening. 

\sa `Kinetic::Ref_counted<T>` 

*/
template< typename Traits, typename Visitor >
class Sort {
public:

/// \name Creation 
/// @{

/*!
The basic constructor. 
*/ 
Sort(Traits tr); 

/// @} 

/// \name Types 
/// @{

/*!
The type of the visitor. 
*/ 
typedef unspecified_type Visitor; 

/*!
The traits type. 
*/ 
typedef unspecified_type Traits; 

/*!
The handle used to refer to vertex in the sorted list. Derefernecing this returns a `Key` into the `ActiveObjectsTable`. 
*/ 
typedef unspecified_type Vertex_handle; 

/*!
A reference counted pointer to be used for storing references to the object. 
*/ 
typedef unspecified_type Handle; 

/*!
A reference counted pointer to be used for storing references to the object. 
*/ 
typedef unspecified_type Const_handle; 

/// @} 

/// \name Operations 
/// @{

/*!
Access the visitor. 
*/ 
Visitor& visitor(); 

/*!
Access the traits. 
*/ 
Traits& traits(); 

/*!
Insert the point. 
*/ 
Vertex_handle insert(Point_key k); 

/*!
Erase the point. 
*/ 
void erase(Vertex_handle k); 

/*!
Swap the pair of objects with `vh` as the first element.  The old
solver `s` is used to compute the next root between the two points
being swapped. This method is called by an Event object.*/
void swap(Vertex_handle vh, typename Traits::Kinetic_kernel::Compare_x_1::result_type& s);

/*!
Begin iterating through the ordered `Vertex_handle`s (the iterator is convertible to `Vertex_handle`. 
*/ 
Iterator begin(); 

/*!
End iterating through the ordered `Vertex_handle`s (the iterator is convertible to `Vertex_handle`. 
*/ 
Iterator end(); 

/// @}

}; /* end Kinetic::Sort */
} /* end namespace Kinetic */
} /* end namespace CGAL */
