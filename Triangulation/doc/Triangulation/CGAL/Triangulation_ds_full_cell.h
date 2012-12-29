
namespace CGAL {

/*!
\ingroup PkgTriangulations

The class `Triangulation_ds_full_cell` serves as the default full cell template parameter in the 
class `Triangulation_data_structure<Dimensionality, TriangulationDSVertex, TriangulationDSFullCell>`. 

This class does not provide any geometric capabilities but only combinatorial 
(adjacency) information. Thus, if the `Triangulation_data_structure` is 
used as a parameter of an (embedded) `Triangulation`, then its full cell template 
parameter has to fulfill additional geometric requirements, i.e. it has to be 
a model of the refined concept `TriangulationFullCell`. 

This class can be used directly or can serve as a base to derive other classes 
with some additional attributes tuned for a specific application. 

Parameters 
-------------- 

The first template parameter, `TriangulationDataStructure`, must be a model of the 
`TriangulationDataStructure` concept. 

The second parameter, `TDSFullCellStoragePolicy`, indicates whether or not 
the full cell should additionally store the mirror indices (the indices 
of the mirror vertices). This improves speed a little, but takes 
more space: 

The class template `Triangulation_ds_full_cell` accepts that no second parameter be specified. 
It also accepts the tag `CGAL::Default` as second parameter. Both cases are 
equivalent to setting `TDSFullCellStoragePolicy` to 
`CGAL::TDS_full_cell_default_storage_policy`. 

When the second parameter is specified, its possible ``values'' 
are:<UL> 

<LI>`CGAL::Default`, which is the default value. In that case, the 
policy `CGAL::TDS_full_cell_default_storage_policy` is used. 

<LI>`CGAL::TDS_full_cell_default_storage_policy`. In that case, the mirror 
indices are not stored. 

<LI>`CGAL::TDS_full_cell_mirror_storage_policy`. In that case, the mirror 
indices are stored. 
</UL> 
See the user manual for how to choose the second option. 

\cgalModels ::TriangulationDSFullCell 

Rebind mechanism 
-------------- 

In case of derivation from that class, the nested class 
`Rebind_TDS` need to be provided in the derived class. 

\sa `Triangulation_ds_vertex<TriangulationDataStructure>` 
\sa `Triangulation_data_structure<Dimensionality, TriangulationDSVertex, TriangulationDSFullCell>>` 

*/
template< typename TriangulationDataStructure, typename TDSFullCellStoragePolicy >
class Triangulation_ds_full_cell {
public:

/// \name Validity check 
/// @{

/*! 

The `is_valid` method is only minimally defined in the
`TriangulationDSFullCell` concept, so that we document it more
precisely here, for the model `Triangulation_ds_full_cell`.

\cgalAdvanced Implements the validity checks required by the concept 
`TriangulationDSFullCell`. In addition, it is checked that there is no 
`NULL` handle to vertices in the middle of non-`NULL` ones, that is, 
that the internal memory layout is not corrupted. 
*/ 
bool is_valid(bool verbose=false) const; 

/// @}

}; /* end Triangulation_ds_full_cell */
} /* end namespace CGAL */
