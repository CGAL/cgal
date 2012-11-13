
/*!
\ingroup PkgLinearCellComplexConcepts
\cgalConcept

The concept `CellAttributeWithPoint` is a refinement of the `CellAttribute` 
concept, to represent a cell attribute containing a point. 

\cgalRefines `CellAttribute` 

\cgalHasModel `CGAL::Cell_attribute_with_point<LCC,Info_,Tag,OnMerge,OnSplit>`

\sa `LinearCellComplexItems` 

*/

class CellAttributeWithPoint {
public:

/// \name Types 
/// @{

/*! 
Type of the used point. 
*/ 
typedef Hidden_type Point; 

/// @} 

/// \name Creation 
/// @{

/*! 
Default constructor. 
*/ 
CellAttributeWithPoint(); 

/*! 
Constructor initializing the point of `cawp` by the 
copy contructor `Point(apoint)`. 
*/ 
CellAttributeWithPoint(const Point&apoint); 

/*! 
Constructor initializing the point of `cawp` by the 
copy contructor `Point(apoint)` and initializing the 
information of `cawp` by the 
copy contructor `Info(info)`. 
Defined only if `Info` is different from `void`. 
*/ 
CellAttributeWithPoint(const Point&apoint, const Info& info); 

/// @} 

/// \name Access Member Functions 
/// @{

/*! 
Returns the point of `cawp`. 
*/ 
Point& point(); 

/*! 
Returns the point of `cawp`, when `cawp` is const. 
*/ 
const Point& point() const; 

/// @}

}; /* end CellAttributeWithPoint */

