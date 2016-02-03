 /*!
\ingroup PkgPointSetProcessingConcepts
\cgalConcept

A PLY interpreter describes how to process a PLY input based on its
headers informations. More specifically, it must specify what
properties are expected for each point and where to store them once
read.

\cgalHasModel `CGAL::Ply_interpreter_points_and_normals_3`

*/
class PlyInterpreter
{
public:
  /*!

    Check if expected properties are found in the PLY hreader.

    \param reader Ply_reader to which the interpreter is expected to
    be applied

    \return true if the provided reader match the requirements of the
    interpreter.
   */
  bool is_applicable (Ply_reader& reader);

  /*!
    Assign read properties to appropriate containers or classes.

    \param reader Ply_reader used to read PLY input
  */
  void process_line (Ply_reader& reader);

};
