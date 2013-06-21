
/*!
\ingroup PkgLinearCellComplexConcepts
\cgalConcept

The concept `CellAttributeWithPoint` is a refinement of the `CellAttribute`
concept, to represent a cell attribute containing a point.

\cgalRefines `CellAttribute`

\cgalHasModel \ref CGAL::Cell_attribute_with_point "CGAL::Cell_attribute_with_point<LCC,Info_,Tag,OnMerge,OnSplit>"

\sa `LinearCellComplexItems`

*/

class CellAttributeWithPoint {
public:

/// \name Types
/// @{

/*!
Type of the used point.
*/
typedef unspecified_type Point;

/*!
Type of the information, defined in the `CellAttribute` concept.
*/
  typedef CellAttribute::Info Info;

/// @}

/// \name Creation
/// @{

/*!
Default constructor.
*/
CellAttributeWithPoint();

/*!
Constructor initializing the point of this attribute by the
copy contructor \ref Point "Point"`(apoint)`.
*/
CellAttributeWithPoint(const Point&apoint);

/*!
Constructor initializing the point of this attribute by the
copy contructor \ref Point "Point"`(apoint)` and initializing the
information of this attribute by the
copy contructor \ref Info "Info"`(info)`.
Defined only if `Info` is different from `void`.
*/
  CellAttributeWithPoint(const Point&apoint, const Info& info);

/// @}

/// \name Access Member Functions
/// @{

/*!
Returns the point of this attribute.
*/
Point& point();

/*!
Returns the point of this attribute, when this is const.
*/
const Point& point() const;

/// @}

}; /* end CellAttributeWithPoint */

