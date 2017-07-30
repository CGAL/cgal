
/*!
\ingroup PkgTSMAConcepts
\cgalConcept

The parameterized shape that is fitted.

\cgalHasModel `PlaneProxy`

*/
class Proxy {
public:
  /// Triangle mesh facet descriptor.
  typedef unspecified_type facet_descriptor;

  /// @name Data members
  /// @{

  /// data member to describe the proxy seed.
  facet_descriptor seed;

  /// }
};

