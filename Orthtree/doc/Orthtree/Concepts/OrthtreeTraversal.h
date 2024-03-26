/*!
  \ingroup PkgOrthtreeConcepts
  \cgalConcept

  \brief Requirements for defining a traversal strategy of an orthtree.

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
    \brief returns the first node of the traversal.
   */
  Node_index first_index() const;

  /*!
    \brief returns the next node to be traversed given the previous node `n`.
   */
  Node_index next(Node_index n) const;
};
