
namespace CGAL {

/*!
\ingroup PkgSTLExtensionUtilities



`Compact` is a tag class. It can be used to parameterize a complexity policy
in order to specify a particularly memory compact variant of an algorithm.
For example, passing `Location_policy<Compact>` as parameter to
`Delaunay_triangulation_3` selects a slower point location which saves memory.



\cgalModels{DefaultConstructible,CopyConstructible}

\sa `Location_policy`
\sa `Fast`
\sa `Fast_location`
\sa `Compact_location`


*/

struct Compact {


}; /* end Compact */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSTLExtensionUtilities



`Fast` is a tag class. It can be used to parameterize a complexity policy
in order to specify a particularly fast variant of an algorithm.
For example, passing `Location_policy<Fast>` as parameter to
`Delaunay_triangulation_3` selects a fast point location at the expense of
more memory usage.



\cgalModels{DefaultConstructible,CopyConstructible}

\sa `Location_policy`
\sa `Compact`
\sa `Compact_location`
\sa `Fast_location`


*/

struct Fast {


}; /* end Fast */
} /* end namespace CGAL */
