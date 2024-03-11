/**
 * @file   ui/gl/ptrs.h
 * @author Gernot Walzl
 * @date   2011-12-19
 */

#ifndef UI_GL_PTRS_H
#define UI_GL_PTRS_H

#include "smarter_ptr.h"

namespace ui { namespace gl {

class Camera;
class OpenGLWindow;
class MainOpenGLWindow;
class KeyboardAdapter;
class MouseAdapter;

typedef SHARED_PTR<Camera> CameraSPtr;
typedef WEAK_PTR<Camera> CameraWPtr;
typedef SHARED_PTR<OpenGLWindow> OpenGLWindowSPtr;
typedef WEAK_PTR<OpenGLWindow> OpenGLWindowWPtr;
typedef SHARED_PTR<MainOpenGLWindow> MainOpenGLWindowSPtr;
typedef WEAK_PTR<MainOpenGLWindow> MainOpenGLWindowWPtr;
typedef SHARED_PTR<KeyboardAdapter> KeyboardAdapterSPtr;
typedef WEAK_PTR<KeyboardAdapter> KeyboardAdapterWPtr;
typedef SHARED_PTR<MouseAdapter> MouseAdapterSPtr;
typedef WEAK_PTR<MouseAdapter> MouseAdapterWPtr;

} }

#endif /* UI_GL_PTRS_H */
