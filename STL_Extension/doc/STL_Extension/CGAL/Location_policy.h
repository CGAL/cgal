
namespace CGAL {

/*!
\ingroup PkgStlExtensionUtilities

A typedef to `Location_policy<Compact>`. 

\sa `Compact` 
\sa `Fast` 
\sa `Location_policy` 
\sa `Fast_location` 

*/
  typedef Location_policy<Compact> Compact_location;
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgStlExtensionUtilities

A typedef to `Location_policy<Fast>`. 

\sa `Compact` 
\sa `Fast` 
\sa `Location_policy` 
\sa `Compact_location` 


*/
  typedef Location_policy<Fast> Fast_location;
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgStlExtensionUtilities



`Location_policy` is a policy class which can be used to specify a trade-off 
between memory usage and time complexity for the point location strategy 
used in a data-structure. 
For example, passing `Location_policy<Compact>` as parameter to 
`Delaunay_triangulation_3` selects a slower point location which saves memory. 



\cgalHeading{Parameters}

`Tag` can only be either `Fast` or `Compact` currently. 

\cgalModels `DefaultConstructible`
\cgalModels `CopyConstructible`

\sa `Compact` 
\sa `Fast` 
\sa `Fast_location` 
\sa `Compact_location` 


*/
template< typename Tag >
class Location_policy {
public:


}; /* end Location_policy */
} /* end namespace CGAL */
