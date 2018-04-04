// SPDX-License-Identifier: GPL-3.0+
#ifndef VIEWER_ACTIONS_H
#define VIEWER_ACTIONS_H
namespace  qglviewer {
/*! Defines the different actions that can be associated with a keyboard
shortcut using setShortcut().

See the <a href="../keyboard.html">keyboard page</a> for details. */
enum KeyboardAction {
  DRAW_AXIS,
  DRAW_GRID,
  DISPLAY_FPS,
  ENABLE_TEXT,
  EXIT_VIEWER,
  CAMERA_MODE,
  FULL_SCREEN,
  STEREO,
  ANIMATION,
  HELP,
  EDIT_CAMERA,
  MOVE_CAMERA_LEFT,
  MOVE_CAMERA_RIGHT,
  MOVE_CAMERA_UP,
  MOVE_CAMERA_DOWN,
  INCREASE_FLYSPEED,
  DECREASE_FLYSPEED
};

/*! Defines the different mouse handlers: camera() or manipulatedFrame().

Used by setMouseBinding(), setMouseBinding(Qt::KeyboardModifiers modifiers,
Qt::MouseButtons, ClickAction, bool, int) and setWheelBinding() to define
which handler receives the mouse events. */
enum MouseHandler { CAMERA, FRAME };

/*! Defines the possible actions that can be binded to a mouse click using
setMouseBinding(Qt::KeyboardModifiers, Qt::MouseButtons, ClickAction, bool,
int).

See the <a href="../mouse.html">mouse page</a> for details. */
enum ClickAction {
  NO_CLICK_ACTION,
  ZOOM_ON_PIXEL,
  ZOOM_TO_FIT,
  SELECT,
  RAP_FROM_PIXEL,
  RAP_IS_CENTER,
  CENTER_FRAME,
  CENTER_SCENE,
  SHOW_ENTIRE_SCENE,
  ALIGN_FRAME,
  ALIGN_CAMERA
};

/*! Defines the possible actions that can be binded to a mouse action (a
click, followed by a mouse displacement).

These actions may be binded to the camera() or to the manipulatedFrame() (see
QGLViewer::MouseHandler) using setMouseBinding(). */
enum MouseAction {
  NO_MOUSE_ACTION,
  ROTATE,
  ZOOM,
  TRANSLATE,
  MOVE_FORWARD,
  LOOK_AROUND,
  MOVE_BACKWARD,
  SCREEN_ROTATE,
  ROLL,
  DRIVE,
  SCREEN_TRANSLATE,
  ZOOM_ON_REGION
};

}
#endif // VIEWER_ACTIONS_H
