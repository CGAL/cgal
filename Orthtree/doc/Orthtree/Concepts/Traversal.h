
/*!
  \ingroup PkgOrthtreeConcepts
  \cgalConcept

  \brief a Traversal provides the functions needed to traverse the
  nodes of an Orthtree.

  A traversal is used to iterate on a tree with a user-selected order
  (e.g. Preorder, Postorder).

  \cgalHasModel `CGAL::Orthtrees::Traversal::Preorder`
  \cgalHasModel `CGAL::Orthtrees::Traversal::Postorder`
  \cgalHasModel `CGAL::Orthtrees::Traversal::Leaves`
 */
class Traversal {
public:

  /*!
    \brief returns the first node to iterate to, given the root of the Orthtree.
   */
  template<typename Node>
  Node first (Node root) const;

  /*!
    \brief returns the next node to iterate to, given the previous node.
   */
  template<typename Node>
  Node next(Node n) const;
};
