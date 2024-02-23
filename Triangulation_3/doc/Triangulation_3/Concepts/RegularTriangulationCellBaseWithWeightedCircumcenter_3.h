
/*!
\ingroup PkgTriangulation3Concepts
\cgalConcept

The concept `RegularTriangulationCellBaseWithWeightedCircumcenter_3` refines
the concept `RegularTriangulationCellBase_3` by caching the weighted circumcenter
of the cell.

\note Since `RegularTriangulationCellBaseWithWeightedCircumcenter_3` refines
the concept `RegularTriangulationCellBase_3`, which requests the function
`weighted_circumcenter()`, the caching mechanism should be implemented directly
into `weighted_circumcenter()`: a `Point_3` is computed when `weighted_circumcenter()` is
first called and its value is stored. In the next calls, the cached value
is returned.

Functions that modify the vertices of the cell will also invalidate the weighted
circumcenter by calling `invalidate_weighted_circumcenter_cache()`.

\cgalRefines{RegularTriangulationCellBase_3}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Regular_triangulation_cell_base_with_weighted_circumcenter_3}
\cgalHasModelsEnd

\sa `RegularTriangulationTraits_3`

*/

class RegularTriangulationCellBaseWithWeightedCircumcenter_3 {
public:

  /*!
    Invalidates the circumcenter value stored in the cell.
  */
  void invalidate_weighted_circumcenter_cache();


}; /* end RegularTriangulationCellBaseWithWeightedCircumcenter_3 */

