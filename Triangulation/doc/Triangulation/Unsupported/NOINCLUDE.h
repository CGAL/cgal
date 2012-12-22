CONVUNSUPPORTED:@ ccRefEnum 
   
  
  
  
  
  
    
om

HEADING:Definition 
--------------

The enum `@` is defined by the class `Triangulation` to specify 
in what kind of face a point has been located in a triangulation.

BEGIN:MEMBER /*!
  \ingroup PkgTriangulations
SPLICE DOC HERE

*/
enum Locate_type {ON_VERTEX, IN_FACE, IN_FACET, IN_FULL_CELL,
OUTSIDE_CONVEX_HULL, OUTSIDE_AFFINE_HULL};

END:MEMBER

HEADING:See Also 
--------------

`CGAL::Triangulation`

CONVUNSUPPORTEDEND
CONVINFO Missing include from ../doc_html/Triangulation/Triangulation_ref/Chapter_main.html.tmp.h
CONVUNSUPPORTED:@ ccRefEnum 

om

Definition 
--------------

The enum `@` is defined by the class `Triangulation` to specify 
in what kind of face a point has been located in a triangulation.

See Also 
--------------

`CGAL::Triangulation`

CONVUNSUPPORTEDEND
CONVINFO Missing include from ../doc_html/Triangulation/Triangulation_ref/Triangulation_face-TriangulationDataStructure-.h

namespace CGAL {

/*!
\ingroup PkgTriangulations

A `Triangulation_face` is a model of the concept `TriangulationDSFace`. 

Parameters 
-------------- 

Parameter `TriangulationDataStructure` must be a model of the concept 
`TriangulationDataStructure`. 
Actually, `Triangulation_face` needs only that this concept defines the types 
`Full_cell_handle`, 
`Vertex_handle`, and 
`Maximal_dimension`. 

\models ::TriangulationDSFace 

\sa `TriangulationDSFace` 
\sa `TriangulationDataStructure` 

*/
template< typename TriangulationDataStructure >
class Triangulation_face {
public:

/// @}

}; /* end Triangulation_face */
} /* end namespace CGAL */
