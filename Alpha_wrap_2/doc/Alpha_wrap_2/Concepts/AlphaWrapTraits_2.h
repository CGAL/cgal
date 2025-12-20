/*!
\ingroup PkgAlphaWrap2Concepts
\cgalConcept

The concept `AlphaWrapTraits_2` defines the requirements for the geometric traits class
of an alpha wrap oracle.

\cgalRefines{AABBRayIntersectionGeomTraits_2, DelaunayTriangulationTraits_2, PolygonTraits_2}

\cgalHasModelsBegin
\cgalHasModelsBare{Any \cgal %kernel is a model of this traits concept}
\cgalHasModelsEnd
*/

class AlphaWrapTraits_2
{
public:
  /*!
  The field type, must be a model of `FieldNumberType` and `FromDoubleConstructible`
  */
  typedef unspecified_type FT;

  /*!
  A predicate object that must provide the following function operator:

  `bool operator()(Segment_2 s)`,

  which returns `true` iff the segment is degenerate.
  */
  typedef unspecified_type Is_degenerate_2;

  // ===

  /*!
  returns the `Is_degenerate_2` predicate.
  */
  Is_degenerate_2 is_degenerate_2_object();
};
