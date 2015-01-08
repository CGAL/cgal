namespace CGAL {

/*!
\ingroup SplitterClasses

Implements the <I>fair</I> splitting rule. 
This splitting rule is a compromise between the median of rectangle 
splitting rule and the midpoint of rectangle splitting rule. This 
splitting rule maintains an upper bound on the maximal allowed ratio 
of the longest and shortest side of a rectangle (the value of this 
upper bound is set in the constructor of the fair splitting 
rule). Among the splits that satisfy this bound, it selects the one in 
which the points have the largest spread. It then splits the points 
in the most even manner possible, subject to maintaining the bound on 
the ratio of the resulting rectangles. 

\cgalHeading{Parameters}

\tparam Traits must be a model of 
the concept `SearchTraits`, for example `CGAL::Search_traits_2`. 

\tparam SpatialSeparator must be a model of the concept `SpatialSeparator`. 
It has as default value the type `Plane_separator<Traits::FT>`. 

\cgalModels `Splitter`

\sa `Splitter` 
\sa `SpatialSeparator` 

*/
template< typename Traits, typename SpatialSeparator >
class Fair {
public:

/// \name Types 
/// @{

/*!
Number type. 
*/ 
typedef Traits::FT FT; 

/// @} 

/// \name Creation 
/// @{

/*!
%Default constructor. 
*/ 
Fair(); 

/*!
Constructor. 
*/ 
Fair(unsigned int bucket_size, 
FT aspect_ratio=FT(3)); 

/// @}

}; /* end Fair */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup SplitterClasses

Implements the <I>median of max spread</I> splitting rule. 
The splitting dimension is the dimension of the longest side of the rectangle. 
The splitting value is defined by the median of the coordinates of the data points 
along this dimension. 

\cgalHeading{Parameters}

Expects for the first template argument a model of 
the concept `SearchTraits`, for example 
the type `CGAL::Search_traits_3< Cartesian<double> >`. 

Expects for the second template argument a model of the concept `SpatialSeparator`. It has as default value 
the type, `CGAL::Plane_separator<Traits::FT>`. 

\cgalModels `Splitter`

\sa `Splitter` 
\sa `SpatialSeparator` 

*/
template< typename Traits, typename SpatialSeparator >
class Median_of_max_spread {
public:

/// \name Creation 
/// @{

/*!
%Default constructor. 
*/ 
Median_of_max_spread(); 

/*!
Constructor. 
*/ 
Median_of_max_spread(unsigned int bucket_size); 

/// @}

}; /* end Median_of_max_spread */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup SplitterClasses

Implements the <I>median of rectangle</I> splitting rule. 
The splitting dimension is the dimension of the longest side of the rectangle. 
The splitting value is defined by the median of the coordinates of the data points 
along this dimension. 

\cgalHeading{Parameters}

Expects for the first template argument a model of 
the concept `SearchTraits`, for example 
the type `CGAL::Search_traits_3< Cartesian<double> >`. 

Expects for the second template argument a model of the concept `SpatialSeparator`. It has as default value 
the type, `CGAL::Plane_separator<Traits::FT>`. 

\cgalModels `Splitter`

\sa `Splitter`
\sa `SpatialSeparator` 

*/
template< typename Traits, typename SpatialSeparator >
class Median_of_rectangle {
public:

/// \name Creation 
/// @{

/*!
%Default constructor. 
*/ 
Median_of_rectangle(); 

/*!
Constructor. 
*/ 
Median_of_rectangle(unsigned int bucket_size); 

/// @}

}; /* end Median_of_rectangle */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup SplitterClasses

Implements the <I>midpoint of max spread</I> splitting rule. 
A rectangle is cut through \f$ (\mathrm{Mind}+\mathrm{Maxd})/2\f$ orthogonal 
to the dimension with the maximum point spread \f$ [\mathrm{Mind},\matrm{Maxd}]\f$. 

\cgalHeading{Parameters}

Expects for the first template argument a model of 
the concept `SearchTraits`, for example 
the type `CGAL::Search_traits_3< Cartesian<double> >`. 

Expects for the second template argument a model of the concept `SpatialSeparator`. It has as default value 
the type, `CGAL::Plane_separator<Traits::FT>` 

\cgalModels `Splitter`

\sa `Splitter`
\sa `SpatialSeparator` 

*/
template< typename Traits, typename SpatialSeparator >
class Midpoint_of_max_spread {
public:

/// \name Creation 
/// @{

/*!
%Default constructor. 
*/ 
Midpoint_of_max_spread(); 

/*!
Constructor. 
*/ 
Midpoint_of_max_spread(unsigned int bucket_size); 

/// @}

}; /* end Midpoint_of_max_spread */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup SplitterClasses

Implements the <I>midpoint of rectangle</I> splitting rule. 
A rectangles is cut through its midpoint orthogonal to the longest side. 

\cgalHeading{Parameters}

Expects for the first template argument a model of 
the concept `SearchTraits`, for example 
the type `CGAL::Search_traits_3< Cartesian<double> >`. 

Expects for the second template argument a model of the concept `SpatialSeparator`. It has as default value 
the type, `CGAL::Plane_separator<Traits::FT>` 

\cgalModels `Splitter`

\sa `Splitter` 
\sa `SpatialSeparator` 

*/
template< typename Traits, typename SpatialSeparator >
class Midpoint_of_rectangle {
public:

/// \name Creation 
/// @{

/*!
%Default constructor. 
*/ 
Midpoint_of_rectangle_spread(); 

/*!
Constructor. 
*/ 
Midpoint_of_rectangle(unsigned int bucket_size); 

/// @}

}; /* end Midpoint_of_rectangle */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup SplitterClasses

Implements the <I>sliding fair</I> splitting rule. 
This splitting rule is a compromise between the `Fair` splitting rule 
and the `Sliding_midpoint` rule. Sliding fair-split is based on the 
theory that there are two types of splits that are good: balanced 
splits that produce fat rectangles, and unbalanced splits provided the 
rectangle with fewer points is fat. 

Also, this splitting rule maintains an upper bound on the maximal 
allowed ratio of the longest and shortest side of a rectangle (the 
value of this upper bound is set in the constructor of the fair 
splitting rule). Among the splits that satisfy this bound, it selects 
the one one in which the points have the largest spread. It then 
considers the most extreme cuts that would be allowed by the aspect 
ratio bound. This is done by dividing the longest side of the 
rectangle by the aspect ratio bound. If the median cut lies between 
these extreme cuts, then we use the median cut. If not, then consider 
the extreme cut that is closer to the median. If all the points lie 
to one side of this cut, then we slide the cut until it hits the first 
point. This may violate the aspect ratio bound, but will never 
generate empty cells. 

\cgalHeading{Parameters}

Expects for the first template argument a model of 
the concept `SearchTraits`, 
for example `CGAL::Cartesian_d<double>`. 

Expects for the second template argument a model of the concept `SpatialSeparator`. It has as default value 
the type, `CGAL::Plane_separator<Traits::FT>` 

\cgalModels `Splitter`

\sa `Splitter` 
\sa `SpatialSeparator` 

*/
template< typename Traits, typename SpatialSeparator >
class Sliding_fair {
public:

/// \name Types 
/// @{

/*!
Number type. 
*/ 
typedef Traits::FT FT; 

/// @} 

/// \name Creation 
/// @{

/*!
Constructor. 
*/ 
Sliding_fair(unsigned int bucket_size, 
FT aspect_ratio=FT(3)); 

/// @} 

/// \name Operations 
/// @{

/*!
Returns the maximal ratio between the largest and smallest side 
of a cell allowed for fair splitting. 
*/ 
FT aspect_ratio(); 

/*!
Returns the bucket size of the leaf nodes. 
*/ 
unsigned int bucket_size(); 

/// @}

}; /* end Sliding_fair */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup SplitterClasses

Implements the <I>sliding midpoint</I> splitting rule. 
This is a modification of the `Midpoint_of_rectangle` splitting rule. 
It first attempts to perform a midpoint of rectangle split as 
described above. If data points lie on both sides of the separating 
plane the sliding midpoint rule computes the same separator as 
the midpoint of rectangle rule. If the data points lie only on one 
side it avoids this by sliding the separator, computed by 
the midpoint of rectangle rule, to the nearest data point. 

\cgalHeading{Parameters}

Expects for the first template argument a model of the concept 
`SearchTraits`, for example `CGAL::Cartesian_d<double>`. 

Expects for the second template argument a model of the concept `SpatialSeparator`. It has as default value 
the type, `CGAL::Plane_separator<Traits::FT>`. 

\cgalModels `Splitter`

\sa `Splitter` 
\sa `SpatialSeparator` 

*/
template< typename Traits, typename SpatialSeparator >
class Sliding_midpoint {
public:

/// \name Creation 
/// @{

/*!
Constructor. 
*/ 
Sliding_midpoint(unsigned int bucket_size); 

/// @} 

/// \name Operations 
/// @{

/*!
Returns the bucket size of the leaf nodes. 
*/ 
unsigned int bucket_size(); 

/// @}

}; /* end Sliding_midpoint */
} /* end namespace CGAL */
