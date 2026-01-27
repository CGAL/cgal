
/*!
\ingroup PkgLinearCellComplexConcepts
\cgalConcept

The concept `CellAttributeWithPoint` is a refinement of the `CellAttribute` concept, to represent a cell attribute containing a point.

\cgalRefines{CellAttribute}

\cgalHasModelsBegin
\cgalHasModelsBare{\link CGAL::Cell_attribute_with_point `CGAL::Cell_attribute_with_point<LCC,Info_,Tag,OnMerge,OnSplit>`\endlink}
\cgalHasModelsEnd

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
copy constructor \link Point `Point`\endlink`(apoint)`.
*/
CellAttributeWithPoint(const Point&apoint);

/*!
Constructor initializing the point of this attribute by the
copy constructor \link Point `Point`\endlink`(apoint)` and initializing the
information of this attribute by the
copy constructor \link Info `Info`\endlink`(info)`.
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

