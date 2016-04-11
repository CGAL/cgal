 /*!
\ingroup PkgPointSetProcessingConcepts
\cgalConcept

A PLY interpreter describes how to process a PLY input based on its
headers informations. More specifically, it must specify the
properties expected for each point and where to store them once read.

\cgalHasModel `CGAL::Ply_interpreter_points_and_normals_3`

*/
class PlyInterpreter
{
public:
  /*!

    Check if expected properties are found in the PLY reader. If this
    function returns false, then `read_ply_custom_points()` will not
    read anything from the input stream.

    Note that although the PLY header specifies the order in which the
    tag attributes should appear in the file, the attributes in the
    `Ply_reader` are only identified by their type and tag, not by
    their position in the header.

    \param reader Ply_reader to which the interpreter is expected to
    be applied

    \return true if the provided reader matches the requirements of
    the interpreter, false otherwise.
   */
  bool is_applicable (Ply_reader& reader);

  /*!

    Assign read properties to appropriate containers or classes. Note
    that in a PLY input, one point occupies exactly one line (with all
    properties separated by white spaces): this method aims at
    processing one point along with its potential attributes.

    \param reader Ply_reader used to read PLY input
  */
  void process_line (Ply_reader& reader);

};
