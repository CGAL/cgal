namespace CGAL {

/*!
  \ingroup PkgGeomviewRef

  An object of the class `Geomview_stream` is a stream in which geometric
  objects can be inserted and where geometric objects can be extracted
  from. The constructor starts Geomview either on the local either on
  a remote machine.

\cgalHeading{Implementation}

  The constructor forks a process and establishes two pipes between the
  processes. The forked process is then overlaid with Geomview. The
  file descriptors `stdin` and `stdout` of Geomview are hooked
  on the two pipes.

  All insert operators construct expressions in `gcl`, the Geomview
  command language, which is a subset of \lisp. These expressions
  are sent to Geomview via the pipe. The extract operators notify `interest`
  for a certain kind of events. When such an event happens Geomview
  sends a description of the event in `gcl` and the extract operator has
  to parse this expression.

  In order to implement further insert and extract operators you should
  take a look at the implementation and at the Geomview manual.

*/
class Geomview_stream {
public:

/// \name Creation
/// @{
  /*!
    Introduces a Geomview stream `gs` with a camera that sees the
    bounding box.  The command `geomview` must be in the user's `PATH`.
    If `machine` and `login` are not `nullptr`,
    Geomview is started on the remote machine using `rsh`.
  */
  Geomview_stream(const Bbox_3 &bbox
                  = Bbox_3(0,0,0, 1,1,1),
                  const char *machine = nullptr,
                  const char *login = nullptr);
/// @}


  /*!
    [`begin`;`end`) is an iterator range with value type
    `Triangle_3<R>`.  This method uses the OFF format to draw several triangles
    at once, which is much faster than drawing them one by one.
  */
  template < class InputIterator >
  void draw_triangles(InputIterator begin, InputIterator end);

/// \name Colors
  /// Geomview distinguishes between edge and face colors. The edge color
  /// is at the same time the color of vertices.

/// @{
  /*!
    Makes `c` the color of vertices, edges and faces in subsequent IO
    operations.
  */
  Geomview_stream& operator<<(const Color& c);

  /*!
    Changes the background color. Returns the old value.
  */
  Color set_bg_color(const Color& c);

  /*!
    Changes the vertex color. Returns the old value.
  */
  Color set_vertex_color(const Color& c);

  /*!
    Changes the edge color. Returns the old value.
  */
  Color set_edge_color(const Color& c);
  /*!
    Changes the face color. Returns the old value.
  */
  Color set_face_color(const Color& c);
/// @}

/// \name Miscellaneous
/// @{
  /*!
    Deletes all objects.
  */
  void clear();

  /*!
    Creates a pickplane (useful after a clear).
  */
  void pickplane();
  /*!
    Positions the camera in a way that all objects can be seen.
  */
  void look_recenter();
  /*!
    Returns the line width.
  */
  int get_line_width() const;
  /*!
    Sets the line width to `w`. Returns the previous value.
  */
  int set_line_width(int w);
  /*!
    Returns the radius of vertices.
  */
  double get_vertex_radius() const;
  /*!
    Sets the radius of vertices to `d`. Returns the previous value.
  */
  double set_vertex_radius(double r) const;
  /*!
    Used to obtain unique identifier names passed to Geomview.  On successive
    calls with the same `s` value, it will return a string which is `s`
    appended with the numbers 0, then 1, then 2...  Note that all counters are
    reset when `clear`() is called.
  */
  string get_new_id(string s);
  /*!
    Sets wired mode.  In wired mode, some structures output only there edges,
    not there surfaces.
    Returns the previous value. By default, wired mode is off.
  */
  bool set_wired(bool b);
  /*!
    Returns `true` iff wired mode is on.
  */
  bool get_wired();

/// @}

/// \name Advanced and Developers Features
/// The following functions are helpful if you develop your own insert
/// and extract functions. The following functions allow to pass a string
/// from Geomview and to read data sent back by Geomview.
/// @{
  /*!
    Inserts string `s` into the stream.
  */
  Geomview_stream& operator<<(string s);

  /*!
    Extracts a string `s` from the stream.
    \pre You have to allocate enough memory.
  */
  Geomview_stream& operator>>(char* s);

  /*!
    Inserts `i` into the stream. Puts whitespace around if the
    stream is in ascii mode.
  */
  Geomview_stream& operator<<(int i);

  /*!
    Inserts `i` into the stream. Puts whitespace around if the
    stream is in ascii mode.
  */
  Geomview_stream& operator<<(unsigned int i);

  /*!
    Inserts `i` into the stream. Puts whitespace around if the
    stream is in ascii mode. Currently implemented by converting to int, so it
    can be truncated on 64 bit platforms.
  */
  Geomview_stream& operator<<(long i);

  /*!
    Inserts `i` into the stream. Puts whitespace around if the
    stream is in ascii mode. Currently implemented by converting to unsigned int,
    so it can be truncated on 64 bit platforms.
  */
  Geomview_stream& operator<<(unsigned long i);

  /*!
    Inserts double `d` into the stream. Puts whitespace around if the
    stream is in ascii mode.
  */
  Geomview_stream& operator<<(double d);

  /*!
    Sets tracing on. The data that are sent to `Geomview` are also
    sent to `cerr`.  Returns the previous value. By default tracing is
    off.
  */
  bool set_trace(bool b);

  /*!
    Returns `true` iff tracing is on.
  */
  bool get_trace();

  /*!
    Sets raw mode.  In raw mode, kernel points are output without headers and
    footers, just the coordinates (in binary or ascii mode).  This allows the
    implementation of the stream functions for other objects to re-use the
    code for points internally, by temporary saving the raw mode to true, and
    restoring it after.
    Returns the previous value. By default, raw mode is off.
  */
  bool set_raw(bool b);

  /*!
    Returns `true` iff raw mode is on.
  */
  bool get_raw();

  /*!
    Sets echo mode.  In echo mode, when you select a point in Geomview, the point
    is actually sent back to Geomview.
    Returns the previous value. By default, echo mode is on.
  */
  bool set_echo(bool b);

  /*!
    Returns `true` iff echo mode is on.
  */
  bool get_echo();

  /*!
    Sets whether we are in binary mode.
  */
  bool set_binary_mode(bool b = true);

  /*!
    Sets whether we are in ascii mode.
  */
  bool set_ascii_mode(bool b = true);

  /*!
    Returns `true` iff `gs` is in binary mode.
  */
  bool get_binary_mode();

  /*!
    Returns `true` iff `gs` is in ascii mode.
  */
  bool get_ascii_mode();

/// @}

}; /* end Geomview_stream */

/*!
  Inserts the point `p` into the stream `gs`.
  \relates Geomview_stream
*/
template <class R>
Geomview_stream&
operator<<(Geomview_stream& gs, const Point_2<R>& p);

/// \addtogroup GeomviewOutput Output Operators for CGAL Kernel Classes
/// \ingroup PkgGeomviewRef
/// The following classes of the \cgal kernel have output
/// operators. 2D objects are embedded in the `xy`-plane.
/// @{

/*!
  Inserts the point `p` into the stream `gs`.
  \relates Geomview_stream
*/
template <class R>
Geomview_stream&
operator<<(Geomview_stream& gs, const Point_3<R>& p);

/*!
  Inserts the segment `s` into the stream `gs`.
  \relates Geomview_stream
*/
template <class R>
Geomview_stream&
operator<<(Geomview_stream& gs, const Segment_2<R>& s);

/*!
  Inserts the segment `s` into the stream `gs`.
  \relates Geomview_stream
*/
template <class R>
Geomview_stream&
operator<<(Geomview_stream& gs, const Segment_3<R>& s);

/*!
  Inserts the ray `r` into the stream `gs`.
  \relates Geomview_stream
*/
template <class R>
Geomview_stream&
operator<<(Geomview_stream& gs, const Ray_2<R>& r);

/*!
  Inserts the ray `r` into the stream `gs`.
  \relates Geomview_stream
*/
template <class R>
Geomview_stream&
operator<<(Geomview_stream& gs, const Ray_3<R>& r);

/*!
  Inserts the line `l` into the stream `gs`.
  \relates Geomview_stream
*/
template <class R>
Geomview_stream&
operator<<(Geomview_stream& gs, const Line_2<R>& l);

/*!
  Inserts the line `l` into the stream `gs`.
  \relates Geomview_stream
*/
template <class R>
Geomview_stream&
operator<<(Geomview_stream& gs, const Line_3<R>& l);

/*!
  Inserts the triangle `t` into the stream `gs`.
  \relates Geomview_stream
*/
template <class R>
Geomview_stream&
operator<<(Geomview_stream& gs, const Triangle_2<R>& t);

/*!
  Inserts the triangle `t` into the stream `gs`.
  \relates Geomview_stream
*/
template <class R>
Geomview_stream&
operator<<(Geomview_stream& gs, const Triangle_3<R>& t);

/*!
  Inserts the tetrahedron `t` into the stream `gs`.
  \relates Geomview_stream
*/
template <class R>
Geomview_stream&
operator<<(Geomview_stream& gs, const Tetrahedron_3<R>& t);

/*!
  Inserts the sphere `s` into the stream `gs`.
  \relates Geomview_stream
*/
template <class R>
Geomview_stream&
operator<<(Geomview_stream& gs, const Sphere_3<R>& s);

/*!
  Inserts the bounding box `b` into the stream `gs`.
  \relates Geomview_stream
*/
Geomview_stream&
operator<<(Geomview_stream& gs, const Bbox_2& b);

/*!
  Inserts the bounding box `b` into the stream `gs`.
  \relates Geomview_stream
*/
Geomview_stream&
operator<<(Geomview_stream& gs, const Bbox_3& b);


/// @}

/// \addtogroup GeomviewInput Input Operators for CGAL Kernel Classes
/// \ingroup PkgGeomviewRef
/// An input operator is provided for points. The user has to select
/// a point on the <I>pick plane</I> with the right mouse button. The pick plane
/// can be moved anywhere with the left mouse button, before a point is entered.
/// @{

/*!
  Extracts the point `p` from the stream `gs`. The point is
  echoed by default, and it depends on the stream echo mode status.
  \relates Geomview_stream
*/
template <class R>
Geomview_stream&
operator>>(Geomview_stream& gs, Point_3<R>& p);

/// @}

/// \addtogroup GeomviewOutputClasses Output Operators for CGAL Basic Library Classes
/// \ingroup PkgGeomviewRef
/// Output operators are provided for polyhedral surfaces, as well as for 3D
/// and 2D triangulations. The latter allow to visualize terrrains if the
/// point type isa  3D point.
/// @{

/*!
  Inserts the polyhedron `P` into the stream `gs`.

  Include: `CGAL/IO/Polyhedron_geomview_ostream.h`

  \relates Geomview_stream
*/
template <class Traits, class HDS>
Geomview_stream&
operator<<(Geomview_stream &G, const Polyhedron_3<Traits,HDS> &P);

/*!

  Inserts the 2D triangulation `T` into the stream `gs`.
  The actual output depends on whether the stream is in wired mode or not.
  Also note that in the case of terrains (when `GT::Point_2` is
  `Point_3<R>`), then the 3D terrain is displayed.

  \relates Geomview_stream
*/
template <class GT, class TDS>
Geomview_stream&
operator<<(Geomview_stream &G, const Triangulation_2<GT,TDS> &T);

/*!
  Inserts the 3D triangulation `T` into the stream `gs`.
  The actual output depends on whether the stream is in wired mode or not.

  \relates Geomview_stream
*/
template <class GT, class TDS>
Geomview_stream&
operator<<(Geomview_stream &G, const Triangulation_3<GT,TDS> &T);

/// @}

} /* end namespace CGAL */
