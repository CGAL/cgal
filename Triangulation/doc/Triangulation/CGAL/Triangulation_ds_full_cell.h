
namespace CGAL {

/*!
\ingroup PkgTriangulationsVertexCellClasses

This class is the default model used for the
full cell of the class `Triangulation_data_structure`.

This class does not provide any geometric capabilities but only combinatorial 
(adjacency) information. Thus, if the `Triangulation_data_structure` is 
used as a parameter of an (embedded) `Triangulation`, then its full cell template 
parameter has to fulfill additional geometric requirements, i.e. it has to be 
a model of the refined concept `TriangulationFullCell`. 

This class can be used directly or can serve as a base to derive other classes 
with some additional attributes tuned for a specific application. 


\tparam TriangulationDataStructure must be a model of the 
`TriangulationDataStructure` concept. 

\tparam TriangulationDSFullCellStoragePolicy indicates whether or not 
the full cell should additionally store the mirror indices (the indices 
of the mirror vertices). This improves speed a little, but takes 
more space.<br>
The class template `Triangulation_ds_full_cell` accepts that no second parameter be specified. 
It also accepts the tag `CGAL::Default` as second parameter. Both cases are 
equivalent to setting `TriangulationDSFullCellStoragePolicy` to 
`CGAL::TDS_full_cell_default_storage_policy`. <br>

When the second parameter is specified, its possible values 
are:<UL> 
<LI>`CGAL::Default`, which is the default value. In that case, the 
policy `CGAL::TDS_full_cell_default_storage_policy` is used (i.e.\ the mirror 
indices are not stored).
<LI>`CGAL::TDS_full_cell_default_storage_policy`. In that case, the mirror 
indices are not stored. 
<LI>`CGAL::TDS_full_cell_mirror_storage_policy`. In that case, the mirror 
indices are stored. 
</UL> 
See the user manual for how to choose the second option. 

\cgalModels `TriangulationDSFullCell`

Rebind mechanism 
-------------- 

In case of derivation from that class, the nested class 
`Rebind_TDS` need to be provided in the derived class. 

\sa `Triangulation_ds_vertex<TriangulationDataStructure>` 
\sa `Triangulation_data_structure<Dimensionality, TriangulationDSVertex, TriangulationDSFullCell>>` 

*/
template< typename TriangulationDataStructure, typename TriangulationDSFullCellStoragePolicy >
class Triangulation_ds_full_cell {
public:

/// \name Validity Check 
/// @{

/*! 
\cgalAdvancedBegin
Implements the validity checks required by the concept 
`TriangulationDSFullCell`.
\cgalAdvancedEnd
*/ 
bool is_valid(bool verbose=false) const; 

/// @}

}; /* end Triangulation_ds_full_cell */
} /* end namespace CGAL */
