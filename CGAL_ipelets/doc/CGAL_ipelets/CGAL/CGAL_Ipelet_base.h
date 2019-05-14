/*!
\ingroup PkgCGALIpelets
The registration of a new ipelet can be done using the macro command `CGAL_IPELET`.
Taking as a parameter the name of the class defining the new ipelet, that macro must be placed in the source
file after the class definition. 
*/
#define CGAL_IPELET(T)

namespace CGAL {

/*!
\ingroup PkgCGALIpelets

`Ipelet_base` is an abstract base class for defining an ipelet. 
The only function that needs to be defined in a derived class is 
`protected_run(int i)` that contains the code of an ipelet. Note that 
the name of the function suggests that the ipelet may throw exceptions that 
will be caught by a function of the class `Ipelet_base` avoiding Ipe to crash. 

\tparam Kernel determines the kernel that will be used
in the ipelet. This parameter must be set according to the algorithm to be used. 
The integer `nbf` indicates the number of functions defined by the ipelet. 


*/
template< typename Kernel, typename int nbf >
class Ipelet_base {
public:

/// \name Types 
/// The set of `read objects` is defined as the set of objects of type
/// `Ipelet_base::Point_2`, `Ipelet_base::Segment_2`, `Ipelet_base::Polygon_2`, `Ipelet_base::Circular_arc_2` or
/// `Ipelet_base::Circle_2`. The set of `drawable objects` is defined as the super set
/// of the set of `read objects` also including objects of type `Ipelet_base::Line_2`,
/// `Ipelet_base::Ray_2`, `Ipelet_base::Triangle_2` and `Ipelet_base::Iso_rectangle_2`.
/// @{

/*!
The kernel used to define internal types. 
*/ 
typedef Kernel Kernel; 

/*!
The number type of the coordinates. 
*/ 
typedef Kernel::FT FT; 

/*!
The point type. 
*/ 
typedef Kernel::Point_2 Point_2; 

/*!
The circle type. 
*/ 
typedef Kernel::Circle_2 Circle_2; 

/*!

The circular arc type. The `CGAL::Sign`, equals either to `CGAL::COUNTERCLOCKWISE` or `CGAL::CLOCKWISE`, indicates 
if the arc is the set of points on the circle from the first point to the second point turning clockwise or counterclockwise. 
*/ 
typedef CGAL::tuple<Circle_2,Point_2,Point_2,CGAL::Sign> Circular_arc_2; 

/*!
The weighted point type. 
*/ 
typedef CGAL::Weighted_point<Point,FT> Weighted_point_2; 

/*!
The segment type. 
*/ 
typedef Kernel::Segment_2 Segment_2; 

/*!
The line type. 
*/ 
typedef Kernel::Line_2 Line_2; 

/*!
The ray type. 
*/ 
typedef Kernel::Ray_2 Ray_2; 

/*!
The triangle type. 
*/ 
typedef Kernel::Triangle_2 Triangle_2; 

/*!
The polygon type. 
*/ 
typedef CGAL::Polygon_2<Kernel> Polygon_2; 

/*!
A type to represent bounding boxes. 
*/ 
typedef Kernel::Iso_rectangle_2 Iso_rectangle_2; 

/*!
Class type providing operators to extract points from segments and polygons. 
*/ 
template<class OutputIterator> Point_grabber; 

/*!
Class type providing operators to extract segments from polygons. 
*/ 
template<class OutputIterator> Segment_grabber; 

/// @} 

/// \name Creation 
/// @{

/*!

initializes an ipelet. 
The string `name` is the name given to the ipelet in the Ipe menu. 
The string array 
`fct_names` contains the names of these functions as they will appear in the sub-menu of the ipelet. 
The string array `help_msg` of size one or `nbf-1` contains a help message for each function of the ipelet. 
These help messages can be printed using the member function `show_help()`. 
This function expects that the last function defined in the ipelet is dedicated to print the help message. 
It is advised that the last function defined is dedicated to show a help message. 

*/ 
CGAL_ipelet(const std::string fct_names[],const std::string help_msg[],const std::string name); 

/// @} 

/// \name Access Functions 
/// @{

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
returns a pointer to an Ipe object representing the drawing page. 
Refer to the Ipe library documentation for more details. 
\cgalAdvancedEnd
*/ 
IpePage* get_ipe_page(); 

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
returns a pointer to an Ipe object providing services to Ipelets. 
Refer to the Ipe library documentation for more details. 
\cgalAdvancedEnd
*/ 
IpeletHelper* get_ipelet_helper(); 

/// @} 

/// \name Operations 
/// @{

/*!
Prints in Ipe a pop-up help message constructed from the string array of the constructor. 
This function expects that the last function defined in the ipelet is dedicated to print the help message. 
When the boolean `one_per_func` is `true`, one help message per function is printed 
(except for the help function itself) using the `nbf-1` strings in the array given to the constructor. Otherwise, 
only one help message is printed using the first string in the array given to the constructor. 

*/ 
void show_help(bool one_per_func=true); 

/*!
Prints the string `msg` as a message in Ipe. 
*/ 
void print_error_message(std::string msg); 

/*!
function called when the user selects the \cgal ipelet. 
*/ 
void protected_run(int i)=0; 

/*!

This function assigns to output iterator `out` all handled objects selected in Ipe. 
Objects read belong to the set of `read objects` (even within groups). 
The output iterator `out` must be able to handle all `read objects` (see `CGAL::Dispatch_or_drop_output_iterator<V,O>`). 
In addition, a bounding box (`Iso_rectangle_2`) of the active objects selected is returned. 
The two Boolean arguments indicate whether the retrieved objects must be deselected and/or removed. 
Note that non-retrieved objects (primitives of Ipe not handled or objects dropped by `out`) 
are automatically deselected, and not deleted. 
If a non-retrieved object is a sub-path or is inside a group, the whole path or group will not be deleted. 

*/ 
template< class V, class O > 
Iso_rectangle_2 read_active_objects(Dispatch_or_drop_output_iterator<V, O> out, 
bool deselect_all=true, 
bool delete_selected_objects=false); 

/*!
returns an output iterator which wraps `it`. `OutputIterator` must be a model 
of the output iterator concept accepting assignments from `Point_2`. 
The returned output iterator will accept assignments from objects of types 
`Polygon_2` or `Segment_2` or `Point_2`, it decomposes them in 
objects of type `Point_2` and assigns them to `it`. 
For more details on the returned output iterator refer to the Boost library 
<A HREF="http://www.boost.org/doc/libs/1_39_0/libs/iterator/doc/function_output_iterator.html">documentation</A>. 

*/ 
template<class OutputIterator> 
boost::function_output_iterator<Point_grabber<OutputIterator> > point_grabber(OutputIterator it); 

/*!
returns an output iterator which wraps `it`. `OutputIterator` must be a model 
of the output iterator concept accepting assignments from `Segment_2`. 
The returned output iterator will accept assignments from objects of types 
`Polygon_2` or `Segment_2`, it decomposes them in objects of type `Segment_2` 
and assigns them to `it`. 
For more details on the returned output iterator refer to the Boost library 
<A HREF="http://www.boost.org/doc/libs/1_39_0/libs/iterator/doc/function_output_iterator.html">documentation</A>. 

*/ 
template<class OutputIterator> 
boost::function_output_iterator<Segment_grabber<OutputIterator> > segment_grabber(OutputIterator it); 

/*!

This function draws in the page of Ipe a given object. 
`T` must be a type inside the set of `drawable objects`. 
When `object` is of type `Line_2` or `Ray_2`, only the part of the object that is inside the page of Ipe (if not empty) is drawn. 
This function is also able to draw a 2D \cgal triangulation as a group of segments. If the boolean `deselect` 
is set to `true`, object drawn is deselected. 

*/ 
template<class T> void draw_in_ipe(const T& object,bool deselect=false); 

/*!

Same as above, except that objects are clipped to `bbox` before been drawn. 

*/ 
template<class T> void draw_in_ipe(const T& object, const Iso_rectangle_2& bbox, bool deselect=false); 

/*!

This function draws in the page of Ipe a set of objects given by an iterator range. 
These objects must be of a type inside the set of `drawable objects`. 
If the boolean `makegrp` is set to `true`, objects drawn define a group in Ipe. If the boolean `deselectall` 
is set to `true`, objects drawn are deselected. 

*/ 
template<class Iterator> void draw_in_ipe(Iterator begin, Iterator end,bool makegrp=true,bool deselectall=false); 

/*!

Same as above, except that objects are clipped to `bbox` before been drawn. 

*/ 
template<class Iterator> void draw_in_ipe(Iterator begin, Iterator end, const Iso_rectangle_2& bbox, bool makegrp=true,bool deselectall=false); 

/*!

This function draws in the page of Ipe a polyline defined by an iterator range of points. 
If the boolean `setclose` is `true`, the polyline drawn is closed. 
If the boolean `deselect` is set to `true`, polyline drawn is deselected. 

*/ 
template<class iterator> void draw_polyline_in_ipe 
(iterator first, iterator last,bool setclose=false,bool deselect=false); 

/*!

This function draws in the page of Ipe the dual of a 2D \cgal triangulation. The edges of the dual are restricted to the interior 
of `bbox`. 
If the boolean `makegrp` is set to `true`, segments drawn define a group in Ipe. If the boolean `deselect` 
is set to `true`, segments drawn are deselected. 

*/ 
template<class Triangulation> void draw_dual_in_ipe(Triangulation T,const Iso_rectangle_2& bbox, 
bool makegrp=true,bool deselect=false); 

/*!

This function draws in the page of Ipe the Voronoi segment skeleton from a triangulation of segments. 
The edges are restricted to the interior of `bbox`. 
If the boolean `makegrp` is set to `true`, segments drawn define a group in Ipe. If the boolean `deselect` 
is set to `true`, segments drawn are deselected. 

*/ 
template<class Triangulation> void draw_skeleton_in_ipe(Triangulation T,const Iso_rectangle_2& bbox, 
bool makegrp=true,bool deselect=false); 

/*!

This function induces the creation of a dialog box requesting a value from the user. 
The string `msg` is printed in this dialog box. 
Ipe lexer tries to convert the user input into an object of type `T` (a simple type such as `int`, `float`, etc). 
If the conversion is possible, the value is returned within a `std::pair`. 
The integer of the pair returned is -1 if the user input is not correct, 0 if user input is empty and 1, otherwise. 
otherwise. 

*/ 
template <class T> std::pair<int,T> request_value_from_user(std::string msg); 

/// @}

}; /* end Ipelet_base */
} /* end namespace CGAL */
