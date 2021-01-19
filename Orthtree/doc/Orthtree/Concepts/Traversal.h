
/*!
  \ingroup PkgOrthtreeConcepts
  \cgalConcept

  \brief a traversal provides the functions needed to traverse the
  nodes of an Orthtree.

  A traversal is used to iterate on a tree with a user-selected order
  (e.g. Preorder, Postorder).

  Template parameters are the template parameters of the `CGAL::Orthtree`,
  please refer to the documentation of this class for more
  information.

  \cgalHasModel `CGAL::Orthtrees::Traversal::Preorder`
  \cgalHasModel `CGAL::Orthtrees::Traversal::Postorder`
  \cgalHasModel `CGAL::Orthtrees::Traversal::Leaves`
 */
template<typename T, typename PR, typename PM>
class Traversal {
public:

  using Node = typename CGAL::Orthtree<T,PR,PM>::Node; ///< The node type

  /*!
    \brief returns the first node to iterate to, given the root of the Orthtree.
   */
  Node first (Node root) const;

  /*!
    \brief returns the next node to iterate to, given the previous node.
   */
  Node next(Node n) const;
};
