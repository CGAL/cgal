namespace CGAL
{

/*!

\ingroup PkgFeatureGraphParameter

The class `Regularization_parameters` describes the parameters for the regularization step.

\sa `CGAL::Regularization_parameters_on_image`
\sa `CGAL::Regularization_parameters_on_surface`
\sa `CGAL::Detect_sharp_features_on_labeled_image`
\sa `CGAL::Detect_sharp_features_on_surface`

*/

class Regularization_parameters
{
public:

  /// \name Types
  /// @{

  /*!
  Numerical type.
  */
  typedef double FT;

  /// @}

  /// \name Access Functions
  /// @{

  /*!
  returns the threshold on the distance between an isthmus line and any line.
  An isthmus line is bounded by a corner that is incident to only the same line.
  In the regularization step, if the maximum distance of a line to another line
  is less than this threshold, then the line is collapsed.
  */
  template <typename Element_type_tag, typename Point_3, typename Index>
  FT get_isthmus_line_distance(const Point_3& point_in_space, const Index& element_index) const;

  /*!
  returns the threshold on the distance between a simple line and any line.
  A simple line is a line that can be removed without splitting the graph of its neighbors.
  In the regularization step, if the maximum distance of a line to another line
  is less than this threshold, then the line is collapsed.
  */
  template <typename Element_type_tag, typename Point_3, typename Index>
  FT get_simple_line_distance(const Point_3& point_in_space, const Index& element_index) const;

  /*!
  returns the threshold on the distance between a corner line and any line.
  A corner line is neither isthmus or simple.
  In the regularization step, if the maximum distance of a line to another line
  is less than this threshold, then the line is collapsed.
  */
  template <typename Element_type_tag, typename Point_3, typename Index>
  FT get_corner_line_distance(const Point_3& point_in_space, const Index& element_index) const;

  /// @}

};

/*!

\ingroup PkgFeatureGraphParameter

The class `Regularization_parameters_on_image` describes the parameters for the regularization step
with default values adapted for image inputs.

\sa `CGAL::Regularization_parameters`
\sa `CGAL::Regularization_parameters_on_surface`
\sa `CGAL::Detect_sharp_features_on_labeled_image`
\sa `CGAL::Detect_sharp_features_on_surface`

*/
class Regularization_parameters_on_image :
public Regularization_parameters
{
private:
  typedef Regularization_parameters Base;

public:

  /// \name Types
  /// @{

  /*!
  Numerical type.
  */
  typedef Base::FT FT;

  /// @}

  /// \name Constructor
  /// @{

  /*!
  constructs the parameters used to regularize the feature graph extracted from an image.

  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below.
            The distances are expressed in terms of the longest voxel edge length.

  \cgalNamedParamsBegin
    \cgalParamSectionBegin{line_distance}
      \cgalParamDescription{the threshold on the distance between lines.
          In the regularization step, if the maximum distance of a line to another line
          is less than this threshold, then the line is collapsed.
          It can be a constant or a functor.
          If it is a functor, it must implement
          <UL>
            <LI> `template <typename Element_type_tag, typename Index>`
            <BR> `FT operator()(const Point_3& point_in_space, const Index& element_index) const`
          </UL>}
      \cgalParamDefault{`FT(4.0)`}
      \cgalParamExtra{This paremeter sets a baseline distance and is overriden
                      by the parameters `parameters::ithmus_line_distance`, `parameters::simple_line_distance` and `parameters::corner_line_distance`
                      if they are set.}
    \cgalParamSectionEnd
    \cgalParamSectionBegin{ithmus_line_distance}
      \cgalParamDescription{the threshold on the distance between an isthmus line and any line.
          An isthmus line is bounded by a corner that is incident to only the same line.
          In the regularization step, if the maximum distance of a line to another line
          is less than this threshold, then the line is collapsed.
          It can be a constant or a functor.
          If it is a functor, it must implement
          <UL>
            <LI> `template <typename Element_type_tag, typename Index>`
            <BR> `FT operator()(const Point_3& point_in_space, const Index& element_index) const`
          </UL>}
      \cgalParamDefault{`FT(-1.0)`}
      \cgalParamExtra{Negative values implies the use of the parameter `parameters::line_distance`.}
    \cgalParamSectionEnd
    \cgalParamSectionBegin{simple_line_distance}
      \cgalParamDescription{the threshold on the distance between a simple line and any line.
          A simple line is a line that can be removed without splitting the graph of its neighbors.
          In the regularization step, if the maximum distance of a line to another line
          is less than this threshold, then the line is collapsed.
          It can be a constant or a functor.
          If it is a functor, it must implement
          <UL>
            <LI> `template <typename Element_type_tag, typename Index>`
            <BR> `FT operator()(const Point_3& point_in_space, const Index& element_index) const`
          </UL>}
      \cgalParamDefault{`FT(-1.0)`}
      \cgalParamExtra{Negative values implies the use of the parameter `parameters::line_distance`.}
    \cgalParamSectionEnd
    \cgalParamSectionBegin{corner_line_distance}
      \cgalParamDescription{the threshold on the distance between a corner line and any line.
          A corner line is neither isthmus or simple.
          In the regularization step, if the maximum distance of a line to another line
          is less than this threshold, then the line is collapsed.
          It can be a constant or a functor.
          If it is a functor, it must implement
          <UL>
            <LI> `template <typename Element_type_tag, typename Index>`
            <BR> `FT operator()(const Point_3& point_in_space, const Index& element_index) const`
          </UL>}
      \cgalParamDefault{`FT(-1.0)`}
      \cgalParamExtra{Negative values implies the use of the parameter `parameters::line_distance`.}
    \cgalParamSectionEnd
  \cgalNamedParamsEnd
  */
  template <typename CGAL_NP_TEMPLATE_PARAMETERS>
  Regularization_parameters_on_image(
    const CGAL_NP_CLASS& np = parameters::default_values());

  /// @}
};


/*!

\ingroup PkgFeatureGraphParameter

The class `Regularization_parameters_on_surface` describes the parameters for the regularization step
with default values adapted for surface inputs.

\sa `CGAL::Regularization_parameters`
\sa `CGAL::Regularization_parameters_on_image`
\sa `CGAL::Detect_sharp_features_on_labeled_image`
\sa `CGAL::Detect_sharp_features_on_surface`

*/
class Regularization_parameters_on_surface :
public Regularization_parameters
{
private:
  typedef Regularization_parameters Base;

public:

  /// \name Types
  /// @{

  /*!
  Numerical type.
  */
  typedef Base::FT FT;

  /// @}

  /// \name Constructor
  /// @{

  /*!
  constructs the parameters used to regularize the feature graph extracted from a surface.

  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below.

  \cgalNamedParamsBegin
    \cgalParamSectionBegin{line_distance}
      \cgalParamDescription{the threshold on the distance between lines.
          In the regularization step, if the maximum distance of a line to another line
          is less than this threshold, then the line is collapsed.
          It can be a constant or a functor.
          If it is a functor, it must implement
          <UL>
            <LI> `template <typename Element_type_tag, typename Index>`
            <BR> `FT operator()(const Point_3& point_in_space, const Index& element_index) const`
          </UL>}
      \cgalParamDefault{`FT(0.0)`}
      \cgalParamExtra{This paremeter sets a baseline distance and is overriden
                      by the parameters `parameters::ithmus_line_distance`, `parameters::simple_line_distance` and `parameters::corner_line_distance`
                      if they are set.}
    \cgalParamSectionEnd
    \cgalParamSectionBegin{ithmus_line_distance}
      \cgalParamDescription{the threshold on the distance between an isthmus line and any line.
          An isthmus line is bounded by a corner that is incident to only the same line.
          In the regularization step, if the maximum distance of a line to another line
          is less than this threshold, then the line is collapsed.
          It can be a constant or a functor.
          If it is a functor, it must implement
          <UL>
            <LI> `template <typename Element_type_tag, typename Index>`
            <BR> `FT operator()(const Point_3& point_in_space, const Index& element_index) const`
          </UL>}
      \cgalParamDefault{`FT(-1.0)`}
      \cgalParamExtra{Negative values implies the use of the parameter `parameters::line_distance`.}
    \cgalParamSectionEnd
    \cgalParamSectionBegin{simple_line_distance}
      \cgalParamDescription{the threshold on the distance between a simple line and any line.
          A simple line is a line that can be removed without splitting the graph of its neighbors.
          In the regularization step, if the maximum distance of a line to another line
          is less than this threshold, then the line is collapsed.
          It can be a constant or a functor.
          If it is a functor, it must implement
          <UL>
            <LI> `template <typename Element_type_tag, typename Index>`
            <BR> `FT operator()(const Point_3& point_in_space, const Index& element_index) const`
          </UL>}
      \cgalParamDefault{`FT(-1.0)`}
      \cgalParamExtra{Negative values implies the use of the parameter `parameters::line_distance`.}
    \cgalParamSectionEnd
    \cgalParamSectionBegin{corner_line_distance}
      \cgalParamDescription{the threshold on the distance between a corner line and any line.
          A corner line is neither isthmus or simple.
          In the regularization step, if the maximum distance of a line to another line
          is less than this threshold, then the line is collapsed.
          It can be a constant or a functor.
          If it is a functor, it must implement
          <UL>
            <LI> `template <typename Element_type_tag, typename Index>`
            <BR> `FT operator()(const Point_3& point_in_space, const Index& element_index) const`
          </UL>}
      \cgalParamDefault{`FT(-1.0)`}
      \cgalParamExtra{Negative values implies the use of the parameter `parameters::line_distance`.}
    \cgalParamSectionEnd
  \cgalNamedParamsEnd
  */
  template <typename CGAL_NP_TEMPLATE_PARAMETERS>
  Regularization_parameters_on_surface(
    const CGAL_NP_CLASS& np = parameters::default_values());

  /// @}
};

} /* namespace CGAL */