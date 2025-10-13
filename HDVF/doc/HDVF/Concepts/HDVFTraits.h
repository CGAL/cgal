/*!
\ingroup PkgHDVFConcepts
\cgalConcept

The concept `HDVFTraits` describes the requirements for geometric traits classes used in this package.

\cgalHasModelsBegin
\cgalHasModelsBare{`CGAL::Homological_discrete_vector_field::HDVF_traits_3<Kernel>`}
\cgalHasModelsEnd
*/
class HDVFTraits
{
public:
    typedef unspecified_type Point ;
    typedef unspecified_type Vector ;
    typedef unspecified_type FT ;
    typedef unspecified_type Bbox ;
};