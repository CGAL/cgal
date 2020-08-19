
/*!
 * \ingroup PkgOctreeConcepts
 * \cgalConcept
 *
 * \brief A Split_criterion is a functor which determines whether an octree node needs to be split during refinement.
 *
 * This can also be done by a function pointer, or a lambda.
 */
struct Split_criterion {


  /// \name Operators
  /// @{

  /*!
   * \brief determines whether a node needs to be split based on a const reference to the node
   *
   * \tparam Node
   * \param n const reference to the node
   * \return whether or not the node needs to be split (true for yes)
   */
  template<class Node>
  bool operator()(const Node &n) const {}

  /// @}
};
