#ifndef DISPLAY_LIST_H
#define DISPLAY_LIST_H

#include <CGAL/gl.h>

struct OpenGL_display_lists_priv;

class OpenGL_display_lists {
  friend struct OpenGL_display_lists_priv;
public:
  OpenGL_display_lists(GLsizei number_of_lists = 1);
  ~OpenGL_display_lists();

  void start(GLsizei list_id = 0, GLenum mode = GL_COMPILE);
  void end(GLsizei list_id = 0);

  void draw(GLsizei list_id = 0);
  GLuint openGL_list_id(GLsizei list_id = 0) const;

  bool begin_draw(GLsizei list_id = 0); /// Calls either draw or start/end,
                                    /// returns true iff the redraw is needed.
  void end_draw(GLsizei list_id = 0);

  bool is_valid(GLsizei list_id = 0) const;
  void invalidate(GLsizei list_id = 0);

protected:
  OpenGL_display_lists_priv* d;
}; // end class OpenGL_display_lists

struct OpenGL_display_lists_priv {
  GLsizei number_of_lists;
  GLuint first_display_list;
  std::vector<bool> validity_vector;
  bool check_lists() {
    if(first_display_list == 0)
    {
      first_display_list = ::glGenLists(number_of_lists);
    }
    return first_display_list != 0; // that would mean an error
  }
}; // end class OpenGL_display_lists_priv

GLuint OpenGL_display_lists::openGL_list_id(GLsizei list_id) const {
  return d->first_display_list + list_id;
}

bool OpenGL_display_lists::is_valid(GLsizei list_id) const {
  return d->validity_vector[list_id];
}

void OpenGL_display_lists::invalidate(GLsizei list_id) {
  d->validity_vector[list_id] = false;
}




OpenGL_display_lists::OpenGL_display_lists(GLsizei number_of_lists)
  : d(new OpenGL_display_lists_priv)
{
  d->number_of_lists = number_of_lists;
  d->first_display_list = 0;
  d->validity_vector.resize(number_of_lists, false);
}

OpenGL_display_lists::~OpenGL_display_lists()
{
  if(d->first_display_list != 0) {
    ::glDeleteLists(d->first_display_list, d->number_of_lists);
  }
  delete d;
}

void OpenGL_display_lists::start(GLsizei list_id, GLenum mode)
{
  if(!d->check_lists()) return;
  ::glNewList(openGL_list_id(list_id), mode);
}

void OpenGL_display_lists::end(GLsizei list_id)
{
  if(!d->check_lists()) return;
  ::glEndList();
  d->validity_vector[list_id] = true;
}

void OpenGL_display_lists::draw(GLsizei list_id)
{
  if(!d->check_lists()) return;
  if(!is_valid(list_id)) return;
  ::glCallList(openGL_list_id(list_id));
}

bool OpenGL_display_lists::begin_draw(GLsizei list_id)
{
  if(is_valid(list_id)) {
    draw(list_id);
    return false;
  }
  else {
    start(list_id, GL_COMPILE);
    return true;
  }
}

void OpenGL_display_lists::end_draw(GLsizei list_id)
{
  if(!is_valid(list_id)) {
    end(list_id);
    draw(list_id);
  }
}

#endif // DISPLAY_LIST_H
