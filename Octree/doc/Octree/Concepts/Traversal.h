#ifndef DOCUMENTATION_TRAVERSAL_H
#define DOCUMENTATION_TRAVERSAL_H

/*!
 * \ingroup PkgOctreeConcepts
 * \cgalConcept
 *
 * \brief Provides the functions needed to traverse the nodes of an octree using a Traversal_iterator.
 *
 * A traversal is used to define a specific walk of the tree (e.g. Preorder, Postorder) iteratively
 * rather than recursively.
 *
 * \todo Link to relevant classes
 * \sa `CGAL::Octree::Traversal_iterator<Value>`
 */
class Traversal {
public:

  /// \name Methods
  /// @{

  /*!
   * \brief uses a reference to the root of the tree to determine the first node of the sequence
   *
   * \todo The template params should be simplified
   * \tparam Point_index
   * \param root a const pointer to the root node
   * \return a const pointer to the first node
   */
  template<class Point_index>
  const Node <Point_index> *first(const Node <Point_index> *root) const {}

  /*!
   * \brief uses a reference to the current node to determine the next node of the sequence
   *
   * \todo The template params should be simplified
   * \tparam Point_index
   * \param n a const pointer to the current node
   * \return a const pointer to the next node
   */
  template<class Point_index>
  const Node <Point_index> *next(const Node <Point_index> *n) const {}

  /// @}

};

#endif //DOCUMENTATION_TRAVERSAL_H
