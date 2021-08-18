/****************************************************************************

 Copyright (c) 2018  GeometryFactory Sarl (France).
 Copyright (C) 2002-2014 Gilles Debunne. All rights reserved.

 This file is part of a fork of the QGLViewer library version 2.7.0.

*****************************************************************************/
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-only

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline

#include <CGAL/license/GraphicsView.h>

#else
#define CGAL_INLINE_FUNCTION
#endif

#include <CGAL/Qt/manipulatedFrame.h>
#include <CGAL/Qt/camera.h>
#include <CGAL/Qt/qglviewer.h>

#include <cstdlib>

#include <QMouseEvent>

namespace CGAL{
namespace qglviewer{


/*! Default constructor.

  The translation is set to (0,0,0), with an identity rotation (0,0,0,1) (see
  Frame constructor for details).

  The different sensitivities are set to their default values (see
  rotationSensitivity(), translationSensitivity(), spinningSensitivity() and
  wheelSensitivity()). */
CGAL_INLINE_FUNCTION
ManipulatedFrame::ManipulatedFrame()
    : action_(NO_MOUSE_ACTION), keepsGrabbingMouse_(false) {
  setRotationSensitivity(1.0);
  setTranslationSensitivity(1.0);
  setSpinningSensitivity(0.3);
  setWheelSensitivity(1.0);
  setZoomSensitivity(1.0);

  isSpinning_ = false;
  previousConstraint_ = nullptr;

  connect(&spinningTimer_, SIGNAL(timeout()), SLOT(spinUpdate()));
}

/*! Equal operator. Calls Frame::operator=() and then copy attributes. */
CGAL_INLINE_FUNCTION
ManipulatedFrame &ManipulatedFrame::operator=(const ManipulatedFrame &mf) {
  Frame::operator=(mf);

  setRotationSensitivity(mf.rotationSensitivity());
  setTranslationSensitivity(mf.translationSensitivity());
  setSpinningSensitivity(mf.spinningSensitivity());
  setWheelSensitivity(mf.wheelSensitivity());
  setZoomSensitivity(mf.zoomSensitivity());

  mouseSpeed_ = 0.0;
  dirIsFixed_ = false;
  keepsGrabbingMouse_ = false;
  action_ = NO_MOUSE_ACTION;

  return *this;
}

/*! Copy constructor. Performs a deep copy of all attributes using operator=().
 */
CGAL_INLINE_FUNCTION
ManipulatedFrame::ManipulatedFrame(const ManipulatedFrame &mf)
    : Frame(mf), MouseGrabber() {
  (*this) = mf;
}

////////////////////////////////////////////////////////////////////////////////

/*! Implementation of the MouseGrabber main method.

The ManipulatedFrame grabsMouse() when the mouse is within a 10 pixels region
around its Camera::projectedCoordinatesOf() position().

See the <a href="../examples/mouseGrabber.html">mouseGrabber example</a> for an
illustration. */
CGAL_INLINE_FUNCTION
void ManipulatedFrame::checkIfGrabsMouse(int x, int y,
                                         const Camera *const camera) {
  const int thresold = 10;
  const Vec proj = camera->projectedCoordinatesOf(position());
  setGrabsMouse(keepsGrabbingMouse_ || ((fabs(x - proj.x) < thresold) &&
                                        (fabs(y - proj.y) < thresold)));
}



////////////////////////////////////////////////////////////////////////////////
//                 M o u s e    h a n d l i n g                               //
////////////////////////////////////////////////////////////////////////////////

/*! Returns \c true when the ManipulatedFrame is being manipulated with the
  mouse.

  Can be used to change the display of the manipulated object during
  manipulation.

  When Camera::frame() of the CGAL::QGLViewer::camera() isManipulated(),
  CGAL::QGLViewer::fastDraw() is used in place of CGAL::QGLViewer::draw() for scene
  rendering. A simplified drawing will then allow for interactive camera
  displacements.  */
CGAL_INLINE_FUNCTION
bool ManipulatedFrame::isManipulated() const {
  return action_ != NO_MOUSE_ACTION;
}

/*! Starts the spinning of the ManipulatedFrame.

This method starts a timer that will call spin() every \p updateInterval
milliseconds. The ManipulatedFrame isSpinning() until you call stopSpinning().
*/
CGAL_INLINE_FUNCTION
void ManipulatedFrame::startSpinning(int updateInterval) {
  isSpinning_ = true;
  spinningTimer_.start(updateInterval);
}

/*! Rotates the ManipulatedFrame by its spinningQuaternion(). Called by a timer
  when the ManipulatedFrame isSpinning(). */
CGAL_INLINE_FUNCTION
void ManipulatedFrame::spin() { rotate(spinningQuaternion()); }

/* spin() and spinUpdate() differ since spin can be used by itself (for instance
   by CGAL::QGLViewer::SCREEN_ROTATE) without a spun emission. Much nicer to use the
   spinningQuaternion() and hence spin() for these incremental updates. Nothing
   special to be done for continuous spinning with this design. */
CGAL_INLINE_FUNCTION
void ManipulatedFrame::spinUpdate() {
  spin();
  Q_EMIT spun();
}

#ifndef DOXYGEN
/*! Protected internal method used to handle mouse events. */
CGAL_INLINE_FUNCTION
void ManipulatedFrame::startAction(int ma, bool withConstraint) {
  action_ = (MouseAction)(ma);

  // #CONNECTION# manipulatedFrame::wheelEvent,
  // manipulatedCameraFrame::wheelEvent and mouseReleaseEvent() restore previous
  // constraint
  if (withConstraint)
    previousConstraint_ = nullptr;
  else {
    previousConstraint_ = constraint();
    setConstraint(nullptr);
  }

  switch (action_) {
  case ROTATE:
  case SCREEN_ROTATE:
    mouseSpeed_ = 0.0;
    stopSpinning();
    break;

  case SCREEN_TRANSLATE:
    dirIsFixed_ = false;
    break;

  default:
    break;
  }
}

/*! Updates mouse speed, measured in pixels/milliseconds. Should be called by
any method which wants to use mouse speed. Currently used to trigger spinning in
mouseReleaseEvent(). */
CGAL_INLINE_FUNCTION
void ManipulatedFrame::computeMouseSpeed(const QMouseEvent *const e) {
  const QPoint delta = (e->pos() - prevPos_);
  const qreal dist = sqrt(qreal(delta.x() * delta.x() + delta.y() * delta.y()));
  delay_ = static_cast<int>(last_move_time.restart());
  if (delay_ == 0)
    // Less than a millisecond: assume delay = 1ms
    mouseSpeed_ = dist;
  else
    mouseSpeed_ = dist / delay_;
}

/*! Return 1 if mouse motion was started horizontally and -1 if it was more
vertical. Returns 0 if this could not be determined yet (perfect diagonal
motion, rare). */
CGAL_INLINE_FUNCTION
int ManipulatedFrame::mouseOriginalDirection(const QMouseEvent *const e) {
  static bool horiz =
      true; // Two simultaneous manipulatedFrame require two mice !

  if (!dirIsFixed_) {
    const QPoint delta = e->pos() - pressPos_;
    dirIsFixed_ = abs(delta.x()) != abs(delta.y());
    horiz = abs(delta.x()) > abs(delta.y());
  }

  if (dirIsFixed_)
    if (horiz)
      return 1;
    else
      return -1;
  else
    return 0;
}

CGAL_INLINE_FUNCTION
qreal ManipulatedFrame::deltaWithPrevPos(QMouseEvent *const event,
                                         Camera *const camera) const {
  qreal dx = qreal(event->x() - prevPos_.x()) / camera->screenWidth();
  qreal dy = qreal(event->y() - prevPos_.y()) / camera->screenHeight();

  qreal value = fabs(dx) > fabs(dy) ? dx : dy;
  return value * zoomSensitivity();
}

CGAL_INLINE_FUNCTION
qreal ManipulatedFrame::wheelDelta(const QWheelEvent *event) const {
  static const qreal WHEEL_SENSITIVITY_COEF = 8E-4;
  return event->angleDelta().y() * wheelSensitivity() * WHEEL_SENSITIVITY_COEF;
}

CGAL_INLINE_FUNCTION
void ManipulatedFrame::zoom(qreal delta, const Camera *const camera) {
  Vec trans(0.0, 0.0, (camera->position() - position()).norm() * delta);

  trans = camera->frame()->orientation().rotate(trans);
  if (referenceFrame())
    trans = referenceFrame()->transformOf(trans);
  translate(trans);
}

#endif // DOXYGEN

/*! Initiates the ManipulatedFrame mouse manipulation.

Overloading of MouseGrabber::mousePressEvent(). See also mouseMoveEvent() and
mouseReleaseEvent().

The mouse behavior depends on which button is pressed. See the <a
href="../mouse.html">CGAL::QGLViewer mouse page</a> for details. */
CGAL_INLINE_FUNCTION
void ManipulatedFrame::mousePressEvent(QMouseEvent *const event,
                                       Camera *const camera) {
  Q_UNUSED(camera);

  if (grabsMouse())
    keepsGrabbingMouse_ = true;

  // #CONNECTION setMouseBinding
  // action_ should no longer possibly be NO_MOUSE_ACTION since this value is
  // not inserted in mouseBinding_
  // if (action_ == NO_MOUSE_ACTION)
  // event->ignore();

  prevPos_ = pressPos_ = event->pos();
}

/*! Modifies the ManipulatedFrame according to the mouse motion.

Actual behavior depends on mouse bindings. See the MouseAction enum
and the <a href="../mouse.html">CGAL::QGLViewer mouse page</a> for details.

The \p camera is used to fit the mouse motion with the display parameters (see
Camera::screenWidth(), Camera::screenHeight(), Camera::fieldOfView()).

Emits the manipulated() signal. */
CGAL_INLINE_FUNCTION
void ManipulatedFrame::mouseMoveEvent(QMouseEvent *const event,
                                      Camera *const camera) {
  switch (action_) {
  case TRANSLATE: {
    const QPoint delta = event->pos() - prevPos_;
    Vec trans(delta.x(), -delta.y(), 0.0);
    // Scale to fit the screen mouse displacement
    switch (camera->type()) {
    case Camera::PERSPECTIVE:
      trans *= 2.0 * tan(camera->fieldOfView() / 2.0) *
               fabs((camera->frame()->coordinatesOf(position())).z) /
               camera->screenHeight();
      break;
    case Camera::ORTHOGRAPHIC: {
      GLdouble w, h;
      camera->getOrthoWidthHeight(w, h);
      trans[0] *= 2.0 * w / camera->screenWidth();
      trans[1] *= 2.0 * h / camera->screenHeight();
      break;
    }
    }
    // Transform to world coordinate system.
    trans =
        camera->frame()->orientation().rotate(translationSensitivity() * trans);
    // And then down to frame
    if (referenceFrame())
      trans = referenceFrame()->transformOf(trans);
    translate(trans);
    break;
  }

  case ZOOM: {
    zoom(deltaWithPrevPos(event, camera), camera);
    break;
  }

  case SCREEN_ROTATE: {
    Vec trans = camera->projectedCoordinatesOf(position());

    const qreal prev_angle =
        atan2(prevPos_.y() - trans[1], prevPos_.x() - trans[0]);
    const qreal angle = atan2(event->y() - trans[1], event->x() - trans[0]);

    const Vec axis =
        transformOf(camera->frame()->inverseTransformOf(Vec(0.0, 0.0, -1.0)));
    Quaternion rot(axis, angle - prev_angle);
    //#CONNECTION# These two methods should go together (spinning detection and
    // activation)
    computeMouseSpeed(event);
    setSpinningQuaternion(rot);
    spin();
    break;
  }

  case SCREEN_TRANSLATE: {
    Vec trans;
    int dir = mouseOriginalDirection(event);
    if (dir == 1)
      trans.setValue(event->x() - prevPos_.x(), 0.0, 0.0);
    else if (dir == -1)
      trans.setValue(0.0, prevPos_.y() - event->y(), 0.0);

    switch (camera->type()) {
    case Camera::PERSPECTIVE:
      trans *= 2.0 * tan(camera->fieldOfView() / 2.0) *
               fabs((camera->frame()->coordinatesOf(position())).z) /
               camera->screenHeight();
      break;
    case Camera::ORTHOGRAPHIC: {
      GLdouble w, h;
      camera->getOrthoWidthHeight(w, h);
      trans[0] *= 2.0 * w / camera->screenWidth();
      trans[1] *= 2.0 * h / camera->screenHeight();
      break;
    }
    }
    // Transform to world coordinate system.
    trans =
        camera->frame()->orientation().rotate(translationSensitivity() * trans);
    // And then down to frame
    if (referenceFrame())
      trans = referenceFrame()->transformOf(trans);

    translate(trans);
    break;
  }

  case ROTATE: {
    Vec trans = camera->projectedCoordinatesOf(position());
    Quaternion rot = deformedBallQuaternion(event->x(), event->y(), trans[0],
                                            trans[1], camera);
    trans = Vec(-rot[0], -rot[1], -rot[2]);
    trans = camera->frame()->orientation().rotate(trans);
    trans = transformOf(trans);
    rot[0] = trans[0];
    rot[1] = trans[1];
    rot[2] = trans[2];
    //#CONNECTION# These two methods should go together (spinning detection and
    // activation)
    computeMouseSpeed(event);
    setSpinningQuaternion(rot);
    spin();
    break;
  }

  case MOVE_FORWARD:
  case MOVE_BACKWARD:
  case LOOK_AROUND:
  case ROLL:
  case DRIVE:
  case ZOOM_ON_REGION:
    // These MouseAction values make no sense for a manipulatedFrame
    break;

  default:
    // Possible when the ManipulatedFrame is a MouseGrabber. This method is then
    // called without startAction because of mouseTracking.
    break;
  }

  if (action_ != NO_MOUSE_ACTION) {
    prevPos_ = event->pos();
    Q_EMIT manipulated();
  }
}

/*! Stops the ManipulatedFrame mouse manipulation.

Overloading of MouseGrabber::mouseReleaseEvent().

If the action was a ROTATE MouseAction, a continuous
spinning is possible if the speed of the mouse cursor is larger than
spinningSensitivity() when the button is released. Press the rotate button again
to stop spinning. See startSpinning() and isSpinning(). */
CGAL_INLINE_FUNCTION
void ManipulatedFrame::mouseReleaseEvent(QMouseEvent *const event,
                                         Camera *const camera) {
  Q_UNUSED(event);
  Q_UNUSED(camera);

  keepsGrabbingMouse_ = false;

  if (previousConstraint_)
    setConstraint(previousConstraint_);

  if (((action_ == ROTATE) ||
       (action_ == SCREEN_ROTATE)) &&
      (mouseSpeed_ >= spinningSensitivity()))
    startSpinning(delay_);

  action_ = NO_MOUSE_ACTION;
}

/*! Overloading of MouseGrabber::mouseDoubleClickEvent().

Left button double click aligns the ManipulatedFrame with the \p camera axis
(see alignWithFrame() and ALIGN_FRAME). Right button projects the
ManipulatedFrame on the \p camera view direction. */
CGAL_INLINE_FUNCTION
void ManipulatedFrame::mouseDoubleClickEvent(QMouseEvent *const event,
                                             Camera *const camera) {
  if (event->modifiers() == ::Qt::NoModifier)
    switch (event->button()) {
    case ::Qt::LeftButton:
      alignWithFrame(camera->frame());
      break;
    case ::Qt::RightButton:
      projectOnLine(camera->position(), camera->viewDirection());
      break;
    default:
      break;
    }
}

/*! Overloading of MouseGrabber::wheelEvent().

Using the wheel is equivalent to a ZOOM MouseAction. See
 CGAL::QGLViewer::setWheelBinding(), setWheelSensitivity(). */
CGAL_INLINE_FUNCTION
void ManipulatedFrame::wheelEvent(QWheelEvent *const event,
                                  Camera *const camera) {
  //#CONNECTION# CGAL::QGLViewer::setWheelBinding
  if (action_ == ZOOM) {
    zoom(wheelDelta(event), camera);
    Q_EMIT manipulated();
  }

  // #CONNECTION# startAction should always be called before
  if (previousConstraint_)
    setConstraint(previousConstraint_);

  action_ = NO_MOUSE_ACTION;
}

////////////////////////////////////////////////////////////////////////////////

/*! Returns "pseudo-distance" from (x,y) to ball of radius size.
\arg for a point inside the ball, it is proportional to the euclidean distance
to the ball \arg for a point outside the ball, it is proportional to the inverse
of this distance (tends to zero) on the ball, the function is continuous. */
static qreal projectOnBall(qreal x, qreal y) {
  // If you change the size value, change angle computation in
  // deformedBallQuaternion().
  const qreal size = 1.0;
  const qreal size2 = size * size;
  const qreal size_limit = size2 * 0.5;

  const qreal d = x * x + y * y;
  return d < size_limit ? sqrt(size2 - d) : size_limit / sqrt(d);
}

#ifndef DOXYGEN
/*! Returns a quaternion computed according to the mouse motion. Mouse positions
are projected on a deformed ball, centered on (\p cx,\p cy). */
Quaternion
CGAL_INLINE_FUNCTION
ManipulatedFrame::deformedBallQuaternion(int x, int y, qreal cx, qreal cy,
                                         const Camera *const camera) {
  // Points on the deformed ball
  qreal px =
      rotationSensitivity() * (prevPos_.x() - cx) / camera->screenWidth();
  qreal py =
      rotationSensitivity() * (cy - prevPos_.y()) / camera->screenHeight();
  qreal dx = rotationSensitivity() * (x - cx) / camera->screenWidth();
  qreal dy = rotationSensitivity() * (cy - y) / camera->screenHeight();

  const Vec p1(px, py, projectOnBall(px, py));
  const Vec p2(dx, dy, projectOnBall(dx, dy));
  // Approximation of rotation angle
  // Should be divided by the projectOnBall size, but it is 1.0
  const Vec axis = cross(p2, p1);
  const qreal angle =
      5.0 *
      asin(sqrt(axis.squaredNorm() / p1.squaredNorm() / p2.squaredNorm()));
  return Quaternion(axis, angle);
}
#endif // DOXYGEN
}}
