// Testing a simple arrangement observer.

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_enums.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

/*
 this program gets a text file that contains operation on an arrangement
 such as inserting and removing curves and points. each operation invokes
 some notifications which compare the expected notification from the input
 file and the actual notification. it is very important to keep the right
 order of notification in order to pass the test.
*/

int ok;
std::ifstream global_input_file;
char one_line[128];
char buff[128];

typedef CGAL::Quotient<CGAL::MP_Float>                Number_type;
typedef CGAL::Cartesian<Number_type>                  Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::X_monotone_curve_2                  Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;

void skip_comments(std::ifstream & is, char* line) {
  while (!is.eof()) {
    is.getline(line, 128);
    if (line[0] != '#') break;
  }
}

void compare_results(std::string str) {
  skip_comments(global_input_file, one_line);
  std::istringstream str_stream(one_line);
  str_stream.getline(buff, 128, ' ');
  if (std::string(buff) != str) {
    std::cout << "Expected " << std::string(buff) << " obtained "
              << str << std::endl;
    ok = -1;
  }
}

// An arrangement observer, used to receive notifications of face splits and
// face mergers.
class Test_observer : public Arrangement_2::Observer {
public:
  using Observer = Arrangement_2::Observer;
  using Base_aos = Arrangement_2::Base_aos;

  using Vertex_handle = Base_aos::Vertex_handle;
  using Halfedge_handle = Base_aos::Halfedge_handle;
  using Face_handle = Base_aos::Face_handle;
  using Ccb_halfedge_circulator = Base_aos::Ccb_halfedge_circulator;
  using Point_2 = Base_aos::Point_2;
  using X_monotone_curve_2 = Base_aos::X_monotone_curve_2;

  Test_observer(Base_aos& arr) : Observer(arr) {}

  /// \name Notification functions on global arrangement operations.
  //@{

  /*! Notification before the arrangement is assigned with the content of
   * another arrangement.
   * \param arr The other arrangement. Notice that the arrangement type is the
   *            type used to instantiate the observer, which is conveniently
   *            defined as `Arrangement_2::Base_aos`.
   */
  virtual void before_assign(const Base_aos& /* arr */) override
  { compare_results("before_assign"); }

  /*!
   * Notification after the arrangement has been assigned with another
   * arrangement.
   */
  virtual void after_assign() override { compare_results("after_assign"); }

  /*! Notification before the arrangement is cleared. */
  virtual void before_clear() override { compare_results("before_clear"); }

  /*!
   * Notification after the arrangement is cleared.
   */
  virtual void after_clear() override { compare_results("after_clear"); }

  /*! Notification before a global operation modifies the arrangement. */
  virtual void before_global_change() override
  { compare_results("before_global_change"); }

  /*! Notification after a global operation is completed. */
  virtual void after_global_change() override
  { compare_results("after_global_change"); }
  //@}

  /// \name Notification functions on observer attachment or detachment.
  //@{

  /*! Notification before the observer is attached to an arrangement.
   * \param arr The arrangement we are about to attach the observer to.
   */
  virtual void before_attach(const Base_aos& /* arr */) override
  { compare_results("before_attach"); }

  /*! Notification after the observer has been attached to an arrangement.
   */
  virtual void after_attach() override { compare_results("after_attach"); }

  /*! Notification before the observer is detached from the arrangement.
   */
  virtual void before_detach() override { compare_results("before_detach"); }

  /*! Notification after the observer has been detached to the arrangement.
   */
  virtual void after_detach() override { compare_results("after_detach"); }
  //@}

  /// \name Notification functions on local changes in the arrangement.
  //@{

  /*! Notification before the creation of a new vertex.
   * \param p The point to be associated with the vertex.
   *          This point cannot lies on the surface boundaries.
   */
  virtual void before_create_vertex(const Point_2& /* p */) override
  { compare_results("before_create_vertex"); }

  /*! Notification after the creation of a new vertex.
   * \param v A handle to the created vertex.
   */
  virtual void after_create_vertex(Vertex_handle /* v */) override
  { compare_results("after_create_vertex"); }

  /*! Notification before the creation of a new boundary vertex.
   * \param cv The curve incident to the surface boundary.
   * \param ind The relevant curve-end.
   * \param bound_x The boundary condition of the vertex in x.
   * \param bound_y The boundary condition of the vertex in y.
   */
  virtual void
  before_create_boundary_vertex(const X_monotone_curve_2& /*cv*/,
                                CGAL::Arr_curve_end /* ind */,
                                CGAL::Arr_parameter_space /* bound_x */,
                                CGAL::Arr_parameter_space /* bound_y */)
    override
  { compare_results("before_create_boundary_vertex"); }

  /*! Notification after the creation of a new vertex at infinity.
   * \param v A handle to the created vertex.
   */
  virtual void after_create_boundary_vertex(Vertex_handle /* v */) override
  { compare_results("after_create_boundary_vertex"); }

  /*! Notification before the creation of a new edge.
   * \param c The x-monotone curve to be associated with the edge.
   * \param v1 A handle to the first end-vertex of the edge.
   * \param v2 A handle to the second end-vertex of the edge.
   */
  virtual void before_create_edge(const X_monotone_curve_2& /* c */,
                                  Vertex_handle /* v1 */,
                                  Vertex_handle /* v2 */) override
  { compare_results("before_create_edge"); }

  /*! Notification after the creation of a new edge.
   * \param e A handle to one of the twin halfedges that were created.
   */
  virtual void after_create_edge(Halfedge_handle /* e */) override
  { compare_results("after_create_edge"); }

  /*! Notification before the modification of an existing vertex.
   * \param v A handle to the vertex to be updated.
   * \param p The point to be associated with the vertex.
   */
  virtual void before_modify_vertex(Vertex_handle /* v */,
                                    const Point_2& /* p */) override
  { compare_results("before_modify_vertex"); }

  /*! Notification after a vertex was modified.
   * \param v A handle to the updated vertex.
   */
  virtual void after_modify_vertex(Vertex_handle /* v */) override
  { compare_results("after_modify_vertex"); }

  /*! Notification before the modification of an existing edge.
   * \param e A handle to one of the twin halfedges to be updated.
   * \param c The x-monotone curve to be associated with the edge.
   */
  virtual void before_modify_edge(Halfedge_handle /* e */,
                                  const X_monotone_curve_2& /* c */) override
  { compare_results("before_modify_edge"); }

  /*! Notification after an edge was modified.
   * \param e A handle to one of the twin halfedges that were updated.
   */
  virtual void after_modify_edge(Halfedge_handle /* e */) override
  { compare_results("after_modify_edge"); }

  /*! Notification before the splitting of an edge into two.
   * \param e A handle to one of the existing halfedges.
   * \param v A vertex representing the split point.
   * \param c1 The x-monotone curve to be associated with the first edge.
   * \param c2 The x-monotone curve to be associated with the second edge.
   */
  virtual void before_split_edge(Halfedge_handle /* e */,
                                 Vertex_handle /* v */,
                                 const X_monotone_curve_2& /* c1 */,
                                 const X_monotone_curve_2& /* c2 */) override
  { compare_results("before_split_edge"); }

  /*! Notification after an edge was split.
   * \param e1 A handle to one of the twin halfedges forming the first edge.
   * \param e2 A handle to one of the twin halfedges forming the second edge.
   */
  virtual void after_split_edge(Halfedge_handle /* e1 */,
                                Halfedge_handle /* e2 */) override
  { compare_results("after_split_edge"); }

  /*! Notification before the splitting of a fictitious edge into two.
   * \param e A handle to one of the existing halfedges.
   * \param v A vertex representing the unbounded split point.
   */
  virtual void before_split_fictitious_edge(Halfedge_handle /* e */,
                                            Vertex_handle /* v */) override
  { compare_results("before_split_fictitious_edge"); }

  /*! Notification after a fictitious edge was split.
   * \param e1 A handle to one of the twin halfedges forming the first edge.
   * \param e2 A handle to one of the twin halfedges forming the second edge.
   */
  virtual void after_split_fictitious_edge(Halfedge_handle /* e1 */,
                                           Halfedge_handle /* e2 */) override
  { compare_results("after_split_fictitious_edge"); }

  /*! Notification before the splitting of a face into two.
   * \param f A handle to the existing face.
   * \param e The new edge whose insertion causes the face to split.
   */
  virtual void before_split_face(Face_handle /* f */,
                                 Halfedge_handle /* e */) override
  { compare_results("before_split_face"); }

  /*! Notification after a face was split.
   * \param f A handle to the face we have just split.
   * \param new_f A handle to the new face that has been created.
   * \param is_hole Whether the new face forms a hole inside f.
   */
  virtual void after_split_face (Face_handle /* f */,
                                 Face_handle /* new_f */,
                                 bool /* is_hole */) override
  { compare_results("after_split_face"); }

  /*! Notification before the splitting of an outer CCB into two.
   * \param f A handle to the face that owns the outer CCB.
   * \param h A circulator representing the component boundary.
   * \param e The new edge whose removal causes the outer CCB to split.
   */
  virtual void before_split_outer_ccb(Face_handle /* f */,
                                      Ccb_halfedge_circulator /* h */,
                                      Halfedge_handle /* e */) override
  { compare_results("before_split_outer_ccb"); }

  /*! Notification after an outer CCB was split.
   * \param f A handle to the face that owns the outer CCBs.
   * \param h1 A circulator representing the boundary of the first component.
   * \param h2 A circulator representing the boundary of the second component.
   */
  virtual void after_split_outer_ccb(Face_handle /* f */,
                                     Ccb_halfedge_circulator /* h1 */,
                                     Ccb_halfedge_circulator /* h2 */) override
  { compare_results("after_split_outer_ccb"); }

  /*! Notification before the splitting of an inner CCB into two.
   * \param f A handle to the face containing the inner CCB.
   * \param h A circulator representing the component boundary.
   * \param e The new edge whose removal causes the inner CCB to split.
   */
  virtual void before_split_inner_ccb(Face_handle /* f */,
                                      Ccb_halfedge_circulator /* h */,
                                      Halfedge_handle /* e */) override
  { compare_results("before_split_inner_ccb"); }

  /*! Notification after an inner CCB was split.
   * \param f A handle to the face containing the inner CCBs.
   * \param h1 A circulator representing the boundary of the first component.
   * \param h2 A circulator representing the boundary of the second component.
   */
  virtual void after_split_inner_ccb(Face_handle /* f */,
                                     Ccb_halfedge_circulator /* h1 */,
                                     Ccb_halfedge_circulator /* h2 */) override
  { compare_results("after_split_inner_ccb"); }

  /*! Notification before the creation of a new outer CCB of a face.
   * \param f A handle to the face that owns the outer CCB.
   * \param e A halfedge along the new outer CCB.
   */
  virtual void before_add_outer_ccb(Face_handle /* f */,
                                    Halfedge_handle /* e */) override
  { compare_results("before_add_outer_ccb"); }

  /*! Notification after an outer CCB was added to a face.
   * \param h A circulator representing the boundary of the new outer CCB.
   */
  virtual void after_add_outer_ccb(Ccb_halfedge_circulator /* h */) override
  { compare_results("after_add_outer_ccb"); }

  /*! Notification before the creation of a new inner CCB inside a face.
   * \param f A handle to the face containing the inner CCB.
   * \param e The new halfedge that forms the new inner CCB.
   */
  virtual void before_add_inner_ccb(Face_handle /* f */,
                                    Halfedge_handle /* e */) override
  { compare_results("before_add_inner_ccb"); }

  /*! Notification after an inner CCB was created inside a face.
   * \param h A circulator representing the boundary of the new inner CCB.
   */
  virtual void after_add_inner_ccb(Ccb_halfedge_circulator /* h */) override
  { compare_results("after_add_inner_ccb"); }

  /*! Notification before the creation of a new isolated vertex inside a face.
   * \param f A handle to the face containing the isolated vertex.
   * \param v The isolated vertex.
   */
  virtual void before_add_isolated_vertex(Face_handle /* f */,
                                          Vertex_handle /* v */) override
  { compare_results("before_add_isolated_vertex"); }

  /*! Notification after an isolated vertex was created inside a face.
   * \param v The isolated vertex.
   */
  virtual void after_add_isolated_vertex(Vertex_handle /* v */) override
  { compare_results("after_add_isolated_vertex"); }

  /*! Notification before the merging of two edges.
   * \param e1 A handle to one of the halfedges forming the first edge.
   * \param e2 A handle to one of the halfedges forming the second edge.
   * \param c The x-monotone curve to be associated with the merged edge.
   */
  virtual void before_merge_edge(Halfedge_handle /* e1 */,
                                 Halfedge_handle /* e2 */,
                                 const X_monotone_curve_2& /* c */) override
  { compare_results("before_merge_edge"); }

  /*! Notification after an edge was merged.
   * \param e A handle to one of the twin halfedges forming the merged edge.
   */
  virtual void after_merge_edge(Halfedge_handle /* e */) override
  { compare_results("after_merge_edge"); }

  /*!
   * Notification before the merging of two fictitious edges.
   * \param e1 A handle to one of the halfedges forming the first edge.
   * \param e2 A handle to one of the halfedges forming the second edge.
   */
  virtual void before_merge_fictitious_edge(Halfedge_handle /* e1 */,
                                            Halfedge_handle /* e2 */) override
  { compare_results("before_merge_fictitious_edge"); }

  /*! Notification after a fictitious edge was merged.
   * \param e A handle to one of the twin halfedges forming the merged edge.
   */
  virtual void after_merge_fictitious_edge (Halfedge_handle /* e */) override
  { compare_results("after_merge_fictitious_edge"); }

  /*! Notification before the merging of two faces.
   * \param f1 A handle to the first face.
   * \param f2 A handle to the second face.
   * \param e The edge whose removal causes the faces to merge.
   */
  virtual void before_merge_face(Face_handle /* f1 */,
                                 Face_handle /* f2 */,
                                 Halfedge_handle /* e */) override
  { compare_results("before_merge_face"); }

  /*! Notification after a face was merged.
   * \param f A handle to the merged face.
   */
  virtual void after_merge_face(Face_handle /* f */) override
  { compare_results("after_merge_face"); }

  /*! Notification before the merging of two outer CCBs.
   * \param f A handle to the face that owns the outer CCBs.
   * \param h1 A circulator representing the boundary of the first component.
   * \param h2 A circulator representing the boundary of the second component.
   * \param e The edge whose insertion or removal causes the CCBs to merge.
   */
  virtual void before_merge_outer_ccb(Face_handle /* f */,
                                      Ccb_halfedge_circulator /* h1 */,
                                      Ccb_halfedge_circulator /* h2 */,
                                      Halfedge_handle /* e */) override
  { compare_results("before_merge_outer_ccb"); }

  /*! Notification after an outer CCB was merged.
   * \param f A handle to the face that owns the outer CCBs.
   * \param h A circulator representing the boundary of the merged component.
   */
  virtual void after_merge_outer_ccb(Face_handle /* f */,
                                     Ccb_halfedge_circulator /* h */) override
  { compare_results("after_merge_outer_ccb"); }

  /*! Notification before the merging of two inner CCBs (holes).
   * \param f A handle to the face that contains the inner CCBs.
   * \param h1 A circulator representing the boundary of the first component.
   * \param h2 A circulator representing the boundary of the second component.
   * \param e The edge whose insertion causes the inner CCBs to merge.
   */
  virtual void before_merge_inner_ccb(Face_handle /* f */,
                                      Ccb_halfedge_circulator /* h1 */,
                                      Ccb_halfedge_circulator /* h2 */,
                                      Halfedge_handle /* e */) override
  { compare_results("before_merge_inner_ccb"); }

  /*! Notification after an inner CCB was merged.
   * \param f A handle to the face that contains the inner CCBs.
   * \param h A circulator representing the boundary of the merged component.
   */
  virtual void after_merge_inner_ccb(Face_handle /* f */,
                                     Ccb_halfedge_circulator /* h */) override
  { compare_results("after_merge_inner_ccb"); }

  /*! Notification before an outer CCB is moved from one face to another.
   * \param from_f A handle to the face that currently owns the outer CCB.
   * \param to_f A handle to the face that should own the outer CCB.
   * \param h A circulator representing the boundary of the component.
   */
  virtual void before_move_outer_ccb(Face_handle /* from_f */,
                                     Face_handle /* to_f */,
                                     Ccb_halfedge_circulator /* h */) override
  { compare_results("before_move_outer_ccb"); }

  /*! Notification after an outer CCB is moved from one face to another.
   * \param h A circulator representing the boundary of the component.
   */
  virtual void after_move_outer_ccb(Ccb_halfedge_circulator /* h */) override
  { compare_results("after_move_outer_ccb"); }


  /*! Notification before an inner CCB is moved from one face to another.
   * \param from_f A handle to the face currently containing the inner CCB.
   * \param to_f A handle to the face that should contain the inner CCB.
   * \param h A circulator representing the boundary of the component.
   */
  virtual void before_move_inner_ccb(Face_handle /* from_f */,
                                     Face_handle /* to_f */,
                                     Ccb_halfedge_circulator /* h */) override
  { compare_results("before_move_inner_ccb"); }

  /*! Notification after an inner CCB is moved from one face to another.
   * \param h A circulator representing the boundary of the component.
   */
  virtual void after_move_inner_ccb(Ccb_halfedge_circulator /* h */) override
  { compare_results("after_move_inner_ccb"); }

  /*! Notification before an isolated vertex is moved from one face to another.
   * \param from_f A handle to the face currently containing the vertex.
   * \param to_f A handle to the face that should contain the vertex.
   * \param v The isolated vertex.
   */
  virtual void before_move_isolated_vertex(Face_handle /* from_f */,
                                           Face_handle /* to_f */,
                                           Vertex_handle /* v */) override
  { compare_results("before_move_isolated_vertex"); }

  /*! Notification after an isolated vertex is moved from one face to another.
   * \param v The isolated vertex.
   */
  virtual void after_move_isolated_vertex(Vertex_handle /* v */) override
  { compare_results("after_move_isolated_vertex"); }

  /*! Notification before the removal of a vertex.
   * \param v A handle to the vertex to be deleted.
   */
  virtual void before_remove_vertex(Vertex_handle /* v */) override
  { compare_results("before_remove_vertex"); }

  /*! Notification after the removal of a vertex.
   */
  virtual void after_remove_vertex() override
  { compare_results("after_remove_vertex"); }

  /*! Notification before the removal of an edge.
   * \param e A handle to one of the twin halfedges to be deleted.
   */
  virtual void before_remove_edge(Halfedge_handle /* e */) override
  { compare_results("before_remove_edge"); }

  /*! Notification after the removal of an edge.
   */
  virtual void after_remove_edge() override
  { compare_results("after_remove_edge"); }

  /*! Notification before the removal of an outer CCB.
   * \param f The face that owns the outer CCB.
   * \param h A circulator representing the boundary of the component.
   */
  virtual void before_remove_outer_ccb(Face_handle /* f */,
                                       Ccb_halfedge_circulator /* h */) override
  { compare_results("before_remove_outer_ccb"); }

  /*! Notification after the removal of an outer CCB.
   * \param f The face that used to own the outer CCB.
   */
  virtual void after_remove_outer_ccb(Face_handle /* f */) override
  { compare_results("after_remove_outer_ccb"); }

  /*! Notification before the removal of an inner CCB.
   * \param f The face containing the inner CCB.
   * \param h A circulator representing the boundary of the component.
   */
  virtual void before_remove_inner_ccb(Face_handle /* f */,
                                       Ccb_halfedge_circulator /* h */) override
  { compare_results("before_remove_inner_ccb"); }

  /*! Notification after the removal of an inner CCB.
   * \param f The face that used to contain the inner CCB.
   */
  virtual void after_remove_inner_ccb(Face_handle /* f */) override
  { compare_results("after_remove_inner_ccb"); }

  //@}

};

int main (int argc, char* argv[]) {
  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << " inputfile" << std::endl;
    CGAL_error();
  }
  else {
    // Construct the arrangement containing one diamond-shaped face.
    Arrangement_2  arr;
    global_input_file.open(argv[1]);
    assert(global_input_file.is_open());
    Test_observer    obs (arr);
    Point_2 p;
    Segment_2 s;
    int i_vh = 0, i_heh = 0;
    std::vector<Arrangement_2::Halfedge_handle> heh_vec;
    std::vector<Arrangement_2::Vertex_handle> vh_vec;
    while (! global_input_file.eof()) {
      skip_comments(global_input_file, one_line);
      std::istringstream str_stream(one_line);
      char c;
      str_stream >> c;
      if (!global_input_file.gcount())
        break;
      if (c=='s') {
        str_stream >> c;
        //read segment
        str_stream >> s;
        if (c=='i') {
          // si means insert intersecting segment
          insert(arr,s);
          std::cout << "intersecting segment insert " << s << std::endl;
        }
        else if (c=='n') {
          // sn means insert non intersecting segment
          heh_vec.push_back(insert_non_intersecting_curve (arr, s));
          std::cout << "non intersecting segment insert " << s << std::endl;
        }
      }
      else if (c=='p') {
        // p means read point
        str_stream >> p ;
        std::cout << "point insert " << p << " index " << vh_vec.size()
                  << std::endl;
        vh_vec.push_back(insert_point(arr,p));
      }
      else if (c=='e') {
        // e means read edge index to be removed
        str_stream >> i_heh ;
        std::cout << "remove edge " << heh_vec[i_heh]->curve() << std::endl;
        remove_edge(arr,heh_vec[i_heh]);
      }
      else if (c=='v') {
        // v means read point index to be removed
        str_stream >> i_vh ;
        std::cout << "remove point " << vh_vec[i_vh]->point() << std::endl;
        remove_vertex(arr,vh_vec[i_vh]);
      }
      else {
        //error
        std::cout << "error, unknowen command" << std::endl;
        return -1;
      }
    }
  }
  return ( ok == 0 ? 0 : -1 );
}
