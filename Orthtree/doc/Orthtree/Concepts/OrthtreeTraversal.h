
/*!
  \ingroup PkgOrthtreeConcepts
  \cgalConcept

  \brief a traversal provides the functions needed to traverse the
  nodes of an orthtree.

  A traversal is used to iterate on a tree with a user-selected order
  (e.g. preorder, postorder).

  \cgalHasModel `CGAL::Orthtrees::Preorder_traversal`
  \cgalHasModel `CGAL::Orthtrees::Postorder_traversal`
  \cgalHasModel `CGAL::Orthtrees::Leaves_traversal`
  \cgalHasModel `CGAL::Orthtrees::Level_traversal`
 */
class OrthtreeTraversal {

public:

  using Node = unspecified_type; ///< A specialization of [CGAL::Orthtree<Traits,PointRange,PointMap>::Node](@ref CGAL::Orthtree::Node)

  /*!
    \brief returns the first node to iterate to, given the root of the Orthtree.
   */
  Node first (Node root) const;

  /*!
    \brief returns the next node to iterate to, given the previous node.
   */
  Node next(Node n) const;
};
