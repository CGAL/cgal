
/*!
  \ingroup PkgOrthtreeConcepts
  \cgalConcept

  \brief a traversal provides the functions needed to traverse the
  nodes of an orthtree.

  A traversal is used to iterate on a tree with a user-selected order
  (e.g. preorder, postorder).

  \cgalHasModelsBegin
  \cgalHasModels{CGAL::Orthtrees::Preorder_traversal}
  \cgalHasModels{CGAL::Orthtrees::Postorder_traversal}
  \cgalHasModels{CGAL::Orthtrees::Leaves_traversal}
  \cgalHasModels{CGAL::Orthtrees::Level_traversal}
  \cgalHasModelsEnd
 */
class OrthtreeTraversal {

public:

  using Node_index = unspecified_type; ///< Index type of the orthtree to be traversed

  /*!
    \brief returns the first node to iterate to, given the root of the Orthtree.
   */
  Node_index first_index() const;

  /*!
    \brief returns the next node to iterate to, given the previous node.
   */
  Node_index next(Node_index n) const;
};
