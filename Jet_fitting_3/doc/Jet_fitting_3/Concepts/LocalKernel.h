
/*!
\ingroup PkgJet_fitting_3Concepts
\cgalConcept

The concept `LocalKernel` describes the set of requirements to be 
fulfilled by any class used to instantiate the second template 
parameter of the class 
`CGAL::Monge_via_jet_fitting<DataKernel,LocalKernel,SvdTraits>`. 

This concept provides the geometric primitives used for the 
computations in the class 
`CGAL::Monge_via_jet_fitting`. 

\cgalHeading{Requirements}

In the class `CGAL::Monge_via_jet_fitting` the scalar type, 
`LocalKernel::FT`, must be the same as that of the `SvdTraits` 
concept : `SvdTraits::FT`. 

The type `LocalKernel::FT` is a model of the FieldWithSqrt concept. 

\cgalHeading{Operations}

The scalar type `LocalKernel::FT` must be a field type with a 
square root. 

Only constructors (from 3 scalars and copy constructors) and access 
methods to coordinates `x()`, `y()`, `z()` are needed for the point and 
vector types. 

\cgalHasModel `CGAL::Cartesian<FieldNumberType>` 
\cgalHasModel `CGAL::Simple_cartesian<FieldNumberType>` 

\sa `DataKernel`
\sa `SvdTraits`

*/

class LocalKernel {
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

}; /* end LocalKernel */

