
/*!
\ingroup PkgPeriodic2Triangulation2Concepts
\cgalConcept

The concept `Periodic_2Offset_2` describes a two-/dimensional integer vector with some specialized access functions and operations. 

\hasModel CGAL::Periodic_2_offset_2 

\sa `Periodic_2TriangulationTraits_2` 
\sa `Periodic_2DelaunayTriangulationTraits_2` 

*/

class Periodic_2Offset_2 {
public:

/// \name Creation 
/// @{

// TODO(NGHK): Check
/*! 
Default constructor. 
*/ 
Periodic_2Offset_2(); 

// TODO(NGHK): Check
/*! 
Constructs the offset (x,y). 
*/ 
Periodic_2Offset_2(int x, int y); 

/// @} 

/// \name Operations 
/// @{

// TODO(NGHK): Check
/*! 
Return the vector sum of `this` and `o`. 
*/ 
Periodic_2Offset_2 operator+(const Periodic_2Offset_2 & o) 
const; 

// TODO(NGHK): Check
/*! 
Return the vector difference of `this` and `o`. 
*/ 
Periodic_2Offset_2 operator-(const Periodic_2Offset_2 & o) 
const; 

// TODO(NGHK): Check
/*! 
Return the negative vector of `this`. 
*/ 
Periodic_2Offset_2 operator-() const; 

// TODO(NGHK): Check
/*! 
Add `o` to `this` using vector addition. 
*/ 
void operator+=(const Periodic_2Offset_2 & o) 
const; 

// TODO(NGHK): Check
/*! 
Subtract `o` from `this` using vector subtraction. 
*/ 
void operator-=(const Periodic_2Offset_2 & o) 
const; 

// TODO(NGHK): Check
/*! 
Return `true` if `o` and `this` represent the same vector. 
*/ 
bool operator==(const Periodic_2Offset_2 & o) const; 

// TODO(NGHK): Check
/*! 
Return `true` if `o` and `this` do not represent the same 
vector. 
*/ 
bool operator!=(const Periodic_2Offset_2 & o) const; 

// TODO(NGHK): Check
/*! 
Compare `this` and `o` lexicographically. 
*/ 
bool operator<(const Periodic_2Offset_2 & o) const; 

/// @} 

/// \name Access Functions 
/// @{

// TODO(NGHK): Check
/*! 
Return the \f$ i\f$-th entry of `this`. 
\pre \f$ i\in\{0,1\}\f$ 
*/ 
int operator[](int i); 

// TODO(NGHK): Check
/*! 
Return the \f$ x\f$-entry of `this`. 
*/ 
int x() const; 

// TODO(NGHK): Check
/*! 
Return the \f$ y\f$-entry of `this`. 
*/ 
int y() const; 

// TODO(NGHK): Check
/*! 
Returns `true` if `this` is equal to (0,0). 
*/ 
bool is_null() const; 

// TODO(NGHK): Check
/*! 
Returns `true` if `this` is equal to (0,0). 
*/ 
bool is_zero() const; 

/// @}

}; /* end Periodic_2Offset_2 */


// TODO(NGHK): Check
/*! 
Inputs an offset from `is`. 
\relates Periodic_2Offset_2 
*/ 
istream& operator>>(istream & is, Periodic_2Offset_2 & off); 

// TODO(NGHK): Check
/*! 
Outputs an offset from `os`. 
\relates Periodic_2Offset_2 
*/ 
ostream& operator<<(ostream & os, Periodic_2Offset_2 & off) const; 
