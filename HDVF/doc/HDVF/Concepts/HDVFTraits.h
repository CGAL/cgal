/*!
\ingroup PkgHDVFConcepts
\cgalConcept

The concept `HDVFTraits` describes the requirements for geometric traits classes used in this package.

\cgalHasModelsBegin
\cgalHasModelsBare{`CGAL::Homological_discrete_vector_field::Hdvf_traits_2<K>`}
\cgalHasModelsBare{`CGAL::Homological_discrete_vector_field::Hdvf_traits_3<K>`}
\cgalHasModelsBare{`CGAL::Homological_discrete_vector_field::Hdvf_traits_d<K>`}
\cgalHasModelsEnd

@todo add requirements on Point, Vector, FT, Bbox
*/
class HDVFTraits
{
public:
    /* \brief Type tag encoding the dimension of the data. Must be a `CGAL::Dimension_tag<>` */
    typedef unspecified_type Dimension ;

    /* \brief Geometric kernel type. Must be a model of the `Kernel` concept. */
    typedef unspecified_type Kernel ;
    typedef unspecified_type Point ;
    typedef unspecified_type Vector ;
    typedef typename Kernel::FT FT ;

    /* \brief Type of bounding boxes in the data dimension. */
    typedef unspecified_type Bbox ;
};