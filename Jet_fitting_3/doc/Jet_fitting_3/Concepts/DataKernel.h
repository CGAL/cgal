
/*!
\ingroup PkgJetFitting3Concepts
\cgalConcept

The concept `DataKernel` describes the set of requirements to be
fulfilled by any class used to instantiate first template parameter of
the class
`CGAL::Monge_via_jet_fitting<DataKernel,LocalKernel,SvdTraits>`.

\cgalHeading{Operations}

Only constructors (from 3 scalars and copy constructors) and access
methods to coordinates `x()`, `y()`, `z()` are needed.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Cartesian<FieldNumberType>}
\cgalHasModels{CGAL::Simple_cartesian<FieldNumberType>}
\cgalHasModelsEnd

\sa `LocalKernel`

*/

class DataKernel {
public:

/// \name Types
/// @{

/*!
The scalar type.
*/
typedef unspecified_type FT;

/*!
The point type.
*/
typedef unspecified_type Point_3;

/*!
The vector type.
*/
typedef unspecified_type Vector_3;

/// @}

}; /* end DataKernel */

