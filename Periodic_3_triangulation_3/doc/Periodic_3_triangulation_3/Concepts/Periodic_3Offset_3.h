
/*!
\ingroup PkgPeriodic3Triangulation3Concepts
\cgalConcept

The concept `Periodic_3Offset_3` describes a three-dimensional integer vector
with some specialized access functions and operations.

\cgalHasModel CGAL::Periodic_3_offset_3

\sa `Periodic_3TriangulationTraits_3`
\sa `Periodic_3DelaunayTriangulationTraits_3`
\sa `Periodic_3RegularTriangulationTraits_3`
*/

class Periodic_3Offset_3 {
public:

/// \name Creation
/// @{

/*!
Default constructor.
*/
Periodic_3Offset_3();

/*!
Constructs the offset (x,y,z).
*/
Periodic_3Offset_3(int x, int y, int z);

/// @}

/// \name Operations
/// @{

/*!
Return the vector sum of `o` and `other`.
*/
Periodic_3Offset_3 operator+(const Periodic_3Offset_3 & other)
const;

/*!
Return the vector difference of `o` and `other`.
*/
Periodic_3Offset_3 operator-(const Periodic_3Offset_3 & other)
const;

/*!
Return the negative vector of `o`.
*/
Periodic_3Offset_3 operator-() const;

/*!
Add `other` to `o` using vector addition.
*/
void operator+=(const Periodic_3Offset_3 & other)
const;

/*!
Subtract `other` from `o` using vector subtraction.
*/
void operator-=(const Periodic_3Offset_3 & other)
const;

/*!
Return `true` if `other` and `o` represent the same vector.
*/
bool operator==(const Periodic_3Offset_3 & other) const;

/*!
Return `true` if `other` and `o` do not represent the same
vector.
*/
bool operator!=(const Periodic_3Offset_3 & other) const;

/*!
Compare `o` and `other` lexicographically.
*/
bool operator<(const Periodic_3Offset_3 & other) const;

/// @}

/// \name Access Functions
/// @{

/*!
Return the \f$ i\f$-th entry of `o`.
*/
int operator[](int i);

/*!
Return the \f$ x\f$-entry of `o`.
*/
int x() const;

/*!
Return the \f$ y\f$-entry of `o`.
*/
int y() const;

/*!
Return the \f$ z\f$-entry of `o`.
*/
int z();

/*!
Returns `true` if `o` is equal to (0,0,0).
*/
bool is_null() const;

/*!
Inputs an offset from `is`.
*/
istream& operator>>(istream & is, Periodic_3_offset_3 & off);

/*!
Outputs an offset from `os`.
*/
ostream& operator<<(ostream & os, Periodic_3_offset_3 & off) const;

/// @}

}; /* end Periodic_3Offset_3 */

