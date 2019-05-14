
namespace CGAL {

/*!
\ingroup PkgIOstreams

An object of the class `Color` is a color available 
for drawing operations in many \cgal output streams. 
Each color is defined by a triple of unsigned chars `(r,g,b)` with 
0 \f$\le\f$ r,g,b \f$ \le \f$ 255, the so-called <I>rgb-value</I> of the color. 

\sa `CGAL::Geomview_stream` 

*/

class Color {
public:

/// \name Creation 
/// @{

/*!
creates a color with rgb-value `(0,0,0)`, i.e.\ black. 
*/ 
Color(); 

/*!
creates a color with rgb-value `(red,green,blue)`. 
*/ 
Color(unsigned char red, unsigned char green, unsigned char blue); 

/// @} 

/// \name Operations 
/// @{

/*!
Test for equality: Two colors are equal, iff their 
rgb-values are equal. 
*/ 
bool operator==(const Color &q) const; 

/*!
Test for inequality. 
*/ 
bool operator!=(const Color &q) const; 

/*!
returns the red component of `c`. 
*/ 
unsigned char red() const; 

/*!
returns the green component of `c`. 
*/ 
unsigned char green() const; 

/*!
returns the blue component of `c`. 
*/ 
unsigned char blue() const; 

/// @} 

/// \name Constants 
/// The following constants are predefined:
/// @{

/*!
Black. 
*/ 
const Color BLACK = Color(0, 0, 0); 

/*!
White. 
*/ 
const Color WHITE = Color(255, 255, 255); 

/*!
Red. 
*/ 
const Color RED = Color(255, 0, 0); 

/*!
Green. 
*/ 
const Color GREEN = Color(0, 255, 0); 

/*!
Blue. 
*/ 
const Color BLUE = Color(0, 0, 255); 

/*!
Violet. 
*/ 
const Color VIOLET = Color(255, 0, 255); 

/*!
Orange. 
*/ 
const Color ORANGE = Color(255, 170, 0); 

/// @}

}; /* end Color */
} /* end namespace CGAL */
