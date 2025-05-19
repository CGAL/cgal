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

#include <CGAL/Qt/frame.h>
#include <CGAL/Qt/constraint.h>
#include <math.h>

namespace CGAL{
namespace qglviewer{


/*! Creates a default Frame.

  Its position() is (0,0,0) and it has an identity orientation() Quaternion. The
  referenceFrame() and the constraint() are \c nullptr. */
CGAL_INLINE_FUNCTION
Frame::Frame() : constraint_(nullptr), referenceFrame_(nullptr) {}

/*! Creates a Frame with a position() and an orientation().

 See the Vec and Quaternion documentations for convenient constructors and
 methods.

 The Frame is defined in the world coordinate system (its referenceFrame() is \c
 nullptr). It has a \c nullptr associated constraint(). */
CGAL_INLINE_FUNCTION
Frame::Frame(const Vec &position, const Quaternion &orientation)
    : t_(position), q_(orientation), constraint_(nullptr), referenceFrame_(nullptr) {}

/*! Equal operator.

  The referenceFrame() and constraint() pointers are copied.

  \attention Signal and slot connections are not copied. */
CGAL_INLINE_FUNCTION
Frame &Frame::operator=(const Frame &frame) {
  // Automatic compiler generated version would not emit the modified() signals
  // as is done in setTranslationAndRotation.
  setTranslationAndRotation(frame.translation(), frame.rotation());
  setConstraint(frame.constraint());
  setReferenceFrame(frame.referenceFrame());
  return *this;
}

/*! Copy constructor.

  The translation() and rotation() as well as constraint() and referenceFrame()
  pointers are copied. */
CGAL_INLINE_FUNCTION
Frame::Frame(const Frame &frame) : QObject() { (*this) = frame; }

/////////////////////////////// MATRICES //////////////////////////////////////

/*! Returns the 4x4 OpenGL transformation matrix represented by the Frame.

  This method should be used in conjunction with \c glMultMatrixd() to modify
  the OpenGL modelview matrix from a Frame hierarchy. With this Frame hierarchy:
  \code
  Frame* body     = new Frame();
  Frame* leftArm  = new Frame();
  Frame* rightArm = new Frame();
  leftArm->setReferenceFrame(body);
  rightArm->setReferenceFrame(body);
  \endcode

  The associated OpenGL drawing code should look like:
  \code
  void Viewer::draw()
  {
        glPushMatrix();
        glMultMatrixd(body->matrix());
        drawBody();

        glPushMatrix();
        glMultMatrixd(leftArm->matrix());
        drawArm();
        glPopMatrix();

        glPushMatrix();
        glMultMatrixd(rightArm->matrix());
        drawArm();
        glPopMatrix();

        glPopMatrix();
  }
  \endcode
  Note the use of nested \c glPushMatrix() and \c glPopMatrix() blocks to
  represent the frame hierarchy: \c leftArm and \c rightArm are both correctly
  drawn with respect to the \c body coordinate system.

  This matrix only represents the local Frame transformation (i.e. with respect
  to the referenceFrame()). Use worldMatrix() to get the full Frame
  transformation matrix (i.e. from the world to the Frame coordinate system).
  These two match when the referenceFrame() is \c nullptr.

  The result is only valid until the next call to matrix(), getMatrix(),
  worldMatrix() or getWorldMatrix(). Use it immediately (as above) or use
  getMatrix() instead.

  \attention The OpenGL format of the result is the transpose of the actual
  mathematical European representation (translation is on the last \e line
  instead of the last \e column).

  \note The scaling factor of the 4x4 matrix is 1.0. */
CGAL_INLINE_FUNCTION
const GLdouble *Frame::matrix() const {
  static GLdouble m[4][4];
  getMatrix(m);
  return (const GLdouble *)(m);
}

/*! \c GLdouble[4][4] version of matrix(). See also getWorldMatrix() and
 * matrix(). */
CGAL_INLINE_FUNCTION
void Frame::getMatrix(GLdouble m[4][4]) const {
  q_.getMatrix(m);

  m[3][0] = t_[0];
  m[3][1] = t_[1];
  m[3][2] = t_[2];
}

/*! \c GLdouble[16] version of matrix(). See also getWorldMatrix() and matrix().
 */
CGAL_INLINE_FUNCTION
void Frame::getMatrix(GLdouble m[16]) const {
  q_.getMatrix(m);

  m[12] = t_[0];
  m[13] = t_[1];
  m[14] = t_[2];
}

/*! Returns a Frame representing the inverse of the Frame space transformation.

  The rotation() of the new Frame is the Quaternion::inverse() of the original
  rotation. Its translation() is the negated inverse rotated image of the
  original translation.

  If a Frame is considered as a space rigid transformation (translation and
  rotation), the inverse() Frame performs the inverse transformation.

  Only the local Frame transformation (i.e. defined with respect to the
  referenceFrame()) is inverted. Use worldInverse() for a global inverse.

  The resulting Frame has the same referenceFrame() as the Frame and a \c nullptr
  constraint().

  \note The scaling factor of the 4x4 matrix is 1.0. */
CGAL_INLINE_FUNCTION
Frame Frame::inverse() const {
  Frame fr(-(q_.inverseRotate(t_)), q_.inverse());
  fr.setReferenceFrame(referenceFrame());
  return fr;
}

/*! Returns the 4x4 OpenGL transformation matrix represented by the Frame.

  This method should be used in conjunction with \c glMultMatrixd() to modify
  the OpenGL modelview matrix from a Frame:
  \code
  // The modelview here corresponds to the world coordinate system.
  Frame fr(pos, Quaternion(from, to));
  glPushMatrix();
  glMultMatrixd(fr.worldMatrix());
  // draw object in the fr coordinate system.
  glPopMatrix();
  \endcode

  This matrix represents the global Frame transformation: the entire
  referenceFrame() hierarchy is taken into account to define the Frame
  transformation from the world coordinate system. Use matrix() to get the local
  Frame transformation matrix (i.e. defined with respect to the
  referenceFrame()). These two match when the referenceFrame() is \c nullptr.

  The OpenGL format of the result is the transpose of the actual mathematical
  European representation (translation is on the last \e line instead of the
  last \e column).

  \attention The result is only valid until the next call to matrix(),
  getMatrix(), worldMatrix() or getWorldMatrix(). Use it immediately (as above)
  or use getWorldMatrix() instead.

  \note The scaling factor of the 4x4 matrix is 1.0. */
CGAL_INLINE_FUNCTION
const GLdouble *Frame::worldMatrix() const {
  // This test is done for efficiency reasons (creates lots of temp objects
  // otherwise).
  if (referenceFrame()) {
    static Frame fr;
    fr.setTranslation(position());
    fr.setRotation(orientation());
    return fr.matrix();
  } else
    return matrix();
}

/*! qreal[4][4] parameter version of worldMatrix(). See also getMatrix() and
 * matrix(). */
CGAL_INLINE_FUNCTION
void Frame::getWorldMatrix(GLdouble m[4][4]) const {
  const GLdouble *mat = worldMatrix();
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      m[i][j] = mat[i * 4 + j];
}

/*! qreal[16] parameter version of worldMatrix(). See also getMatrix() and
 * matrix(). */
CGAL_INLINE_FUNCTION
void Frame::getWorldMatrix(GLdouble m[16]) const {
  const GLdouble *mat = worldMatrix();
  for (int i = 0; i < 16; ++i)
    m[i] = mat[i];
}

/*! This is an overloaded method provided for convenience. Same as
 * setFromMatrix(). */
CGAL_INLINE_FUNCTION
void Frame::setFromMatrix(const GLdouble m[4][4]) {
  if (fabs(m[3][3]) < 1E-8) {
    qWarning("Frame::setFromMatrix: Null homogeneous coefficient");
    return;
  }

  qreal rot[3][3];
  for (int i = 0; i < 3; ++i) {
    t_[i] = m[3][i] / m[3][3];
    for (int j = 0; j < 3; ++j)
      // Beware of the transposition (OpenGL to European math)
      rot[i][j] = m[j][i] / m[3][3];
  }
  q_.setFromRotationMatrix(rot);
  Q_EMIT modified();
}

/*! Sets the Frame from an OpenGL matrix representation (rotation in the upper
 left 3x3 matrix and translation on the last line).

 Hence, if a code fragment looks like:
 \code
 GLdouble m[16]={...};
 glMultMatrixd(m);
 \endcode
 It is equivalent to write:
 \code
 Frame fr;
 fr.setFromMatrix(m);
 glMultMatrixd(fr.matrix());
 \endcode

 Using this conversion, you can benefit from the powerful Frame transformation
 methods to translate points and vectors to and from the Frame coordinate system
 to any other Frame coordinate system (including the world coordinate system).
 See coordinatesOf() and transformOf().

 Emits the modified() signal. See also matrix(), getMatrix() and
 Quaternion::setFromRotationMatrix().

 \attention A Frame does not contain a scale factor. The possible scaling in \p
 m will not be converted into the Frame by this method. */
CGAL_INLINE_FUNCTION
void Frame::setFromMatrix(const GLdouble m[16]) {
  GLdouble mat[4][4];
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      mat[i][j] = m[i * 4 + j];
  setFromMatrix(mat);
}

//////////////////// SET AND GET LOCAL TRANSLATION AND ROTATION
//////////////////////////////////

/*! Same as setTranslation(), but with \p qreal parameters. */
CGAL_INLINE_FUNCTION
void Frame::setTranslation(qreal x, qreal y, qreal z) {
  setTranslation(Vec(x, y, z));
}

/*! Fill \c x, \c y and \c z with the translation() of the Frame. */
CGAL_INLINE_FUNCTION
void Frame::getTranslation(qreal &x, qreal &y, qreal &z) const {
  const Vec t = translation();
  x = t[0];
  y = t[1];
  z = t[2];
}

/*! Same as setRotation() but with \c qreal Quaternion parameters. */
CGAL_INLINE_FUNCTION
void Frame::setRotation(qreal q0, qreal q1, qreal q2, qreal q3) {
  setRotation(Quaternion(q0, q1, q2, q3));
}

/*! The \p q are set to the rotation() of the Frame.

See Quaternion::Quaternion(qreal, qreal, qreal, qreal) for details on \c q. */
CGAL_INLINE_FUNCTION
void Frame::getRotation(qreal &q0, qreal &q1, qreal &q2, qreal &q3) const {
  const Quaternion q = rotation();
  q0 = q[0];
  q1 = q[1];
  q2 = q[2];
  q3 = q[3];
}

////////////////////////////////////////////////////////////////////////////////

/*! Translates the Frame of \p t (defined in the Frame coordinate system).

  The translation actually applied to the Frame may differ from \p t since it
  can be filtered by the constraint(). Use translate(Vec&) or
  setTranslationWithConstraint() to retrieve the filtered translation value. Use
  setTranslation() to directly translate the Frame without taking the
  constraint() into account.

  See also rotate(const Quaternion&). Emits the modified() signal. */
CGAL_INLINE_FUNCTION
void Frame::translate(const Vec &t) {
  Vec tbis = t;
  translate(tbis);
}

/*! Same as translate(const Vec&) but \p t may be modified to satisfy the
  translation constraint(). Its new value corresponds to the translation that
  has actually been applied to the Frame. */
CGAL_INLINE_FUNCTION
void Frame::translate(Vec &t) {
  if (constraint())
    constraint()->constrainTranslation(t, this);
  t_ += t;
  Q_EMIT modified();
}

/*! Same as translate(const Vec&) but with \c qreal parameters. */
CGAL_INLINE_FUNCTION
void Frame::translate(qreal x, qreal y, qreal z) {
  Vec t(x, y, z);
  translate(t);
}

/*! Same as translate(Vec&) but with \c qreal parameters. */
CGAL_INLINE_FUNCTION
void Frame::translate(qreal &x, qreal &y, qreal &z) {
  Vec t(x, y, z);
  translate(t);
  x = t[0];
  y = t[1];
  z = t[2];
}

/*! Rotates the Frame by \p q (defined in the Frame coordinate system): R = R*q.

  The rotation actually applied to the Frame may differ from \p q since it can
  be filtered by the constraint(). Use rotate(Quaternion&) or
  setRotationWithConstraint() to retrieve the filtered rotation value. Use
  setRotation() to directly rotate the Frame without taking the constraint()
  into account.

  See also translate(const Vec&). Emits the modified() signal. */
CGAL_INLINE_FUNCTION
void Frame::rotate(const Quaternion &q) {
  Quaternion qbis = q;
  rotate(qbis);
}

/*! Same as rotate(const Quaternion&) but \p q may be modified to satisfy the
  rotation constraint(). Its new value corresponds to the rotation that has
  actually been applied to the Frame. */
CGAL_INLINE_FUNCTION
void Frame::rotate(Quaternion &q) {
  if (constraint())
    constraint()->constrainRotation(q, this);
  q_ *= q;
  q_.normalize(); // Prevents numerical drift
  Q_EMIT modified();
}

/*! Same as rotate(Quaternion&) but with \c qreal Quaternion parameters. */
CGAL_INLINE_FUNCTION
void Frame::rotate(qreal &q0, qreal &q1, qreal &q2, qreal &q3) {
  Quaternion q(q0, q1, q2, q3);
  rotate(q);
  q0 = q[0];
  q1 = q[1];
  q2 = q[2];
  q3 = q[3];
}

/*! Same as rotate(const Quaternion&) but with \c qreal Quaternion parameters.
 */
CGAL_INLINE_FUNCTION
void Frame::rotate(qreal q0, qreal q1, qreal q2, qreal q3) {
  Quaternion q(q0, q1, q2, q3);
  rotate(q);
}

/*! Makes the Frame rotate() by \p rotation around \p point.

  \p point is defined in the world coordinate system, while the \p rotation axis
  is defined in the Frame coordinate system.

  If the Frame has a constraint(), \p rotation is first constrained using
  Constraint::constrainRotation(). The translation which results from the
  filtered rotation around \p point is then computed and filtered using
  Constraint::constrainTranslation(). The new \p rotation value corresponds to
  the rotation that has actually been applied to the Frame.

  Emits the modified() signal. */
CGAL_INLINE_FUNCTION
void Frame::rotateAroundPoint(Quaternion &rotation, const Vec &point) {
  if (constraint())
    constraint()->constrainRotation(rotation, this);
  q_ *= rotation;
  q_.normalize(); // Prevents numerical drift
  Vec trans = point +
              Quaternion(inverseTransformOf(rotation.axis()), rotation.angle())
                  .rotate(position() - point) -
              t_;
  if (constraint())
    constraint()->constrainTranslation(trans, this);
  t_ += trans;
  Q_EMIT modified();
}

/*! Same as rotateAroundPoint(), but with a \c const \p rotation Quaternion.
  Note that the actual rotation may differ since it can be filtered by the
  constraint(). */
CGAL_INLINE_FUNCTION
void Frame::rotateAroundPoint(const Quaternion &rotation, const Vec &point) {
  Quaternion rot = rotation;
  rotateAroundPoint(rot, point);
}

//////////////////// SET AND GET WORLD POSITION AND ORIENTATION
//////////////////////////////////

/*! Sets the position() of the Frame, defined in the world coordinate system.
Emits the modified() signal.

Use setTranslation() to define the \e local frame translation (with respect to
the referenceFrame()). The potential constraint() of the Frame is not taken into
account, use setPositionWithConstraint() instead. */
CGAL_INLINE_FUNCTION
void Frame::setPosition(const Vec &position) {
  if (referenceFrame())
    setTranslation(referenceFrame()->coordinatesOf(position));
  else
    setTranslation(position);
}

/*! Same as setPosition(), but with \c qreal parameters. */
CGAL_INLINE_FUNCTION
void Frame::setPosition(qreal x, qreal y, qreal z) {
  setPosition(Vec(x, y, z));
}

/*! Same as successive calls to setPosition() and then setOrientation().

Only one modified() signal is emitted, which is convenient if this signal is
connected to a CGAL::QGLViewer::update() slot. See also setTranslationAndRotation()
and setPositionAndOrientationWithConstraint(). */
CGAL_INLINE_FUNCTION
void Frame::setPositionAndOrientation(const Vec &position,
                                      const Quaternion &orientation) {
  if (referenceFrame()) {
    t_ = referenceFrame()->coordinatesOf(position);
    q_ = referenceFrame()->orientation().inverse() * orientation;
  } else {
    t_ = position;
    q_ = orientation;
  }
  Q_EMIT modified();
}

/*! Same as successive calls to setTranslation() and then setRotation().

Only one modified() signal is emitted, which is convenient if this signal is
connected to a CGAL::QGLViewer::update() slot. See also setPositionAndOrientation()
and setTranslationAndRotationWithConstraint(). */
CGAL_INLINE_FUNCTION
void Frame::setTranslationAndRotation(const Vec &translation,
                                      const Quaternion &rotation) {
  t_ = translation;
  q_ = rotation;
  Q_EMIT modified();
}

/*! \p x, \p y and \p z are set to the position() of the Frame. */
CGAL_INLINE_FUNCTION
void Frame::getPosition(qreal &x, qreal &y, qreal &z) const {
  Vec p = position();
  x = p.x;
  y = p.y;
  z = p.z;
}

/*! Sets the orientation() of the Frame, defined in the world coordinate system.
Emits the modified() signal.

Use setRotation() to define the \e local frame rotation (with respect to the
referenceFrame()). The potential constraint() of the Frame is not taken into
account, use setOrientationWithConstraint() instead. */
CGAL_INLINE_FUNCTION
void Frame::setOrientation(const Quaternion &orientation) {
  if (referenceFrame())
    setRotation(referenceFrame()->orientation().inverse() * orientation);
  else
    setRotation(orientation);
}

/*! Same as setOrientation(), but with \c qreal parameters. */
CGAL_INLINE_FUNCTION
void Frame::setOrientation(qreal q0, qreal q1, qreal q2, qreal q3) {
  setOrientation(Quaternion(q0, q1, q2, q3));
}

/*! Get the current orientation of the frame (same as orientation()).
  Parameters are the orientation Quaternion values.
  See also setOrientation(). */

/*! The \p q are set to the orientation() of the Frame.

See Quaternion::Quaternion(qreal, qreal, qreal, qreal) for details on \c q. */
CGAL_INLINE_FUNCTION
void Frame::getOrientation(qreal &q0, qreal &q1, qreal &q2, qreal &q3) const {
  Quaternion o = orientation();
  q0 = o[0];
  q1 = o[1];
  q2 = o[2];
  q3 = o[3];
}

/*! Returns the position of the Frame, defined in the world coordinate system.
   See also orientation(), setPosition() and translation(). */
CGAL_INLINE_FUNCTION
Vec Frame::position() const {
  if (referenceFrame_)
    return inverseCoordinatesOf(Vec(0.0, 0.0, 0.0));
  else
    return t_;
}

/*! Returns the orientation of the Frame, defined in the world coordinate
  system. See also position(), setOrientation() and rotation(). */
CGAL_INLINE_FUNCTION
Quaternion Frame::orientation() const {
  Quaternion res = rotation();
  const Frame *fr = referenceFrame();
  while (fr != nullptr) {
    res = fr->rotation() * res;
    fr = fr->referenceFrame();
  }
  return res;
}

////////////////////// C o n s t r a i n t   V e r s i o n s
/////////////////////////////

/*! Same as setTranslation(), but \p translation is modified so that the
  potential constraint() of the Frame is satisfied.

  Emits the modified() signal. See also setRotationWithConstraint() and
  setPositionWithConstraint(). */
CGAL_INLINE_FUNCTION
void Frame::setTranslationWithConstraint(Vec &translation) {
  Vec deltaT = translation - this->translation();
  if (constraint())
    constraint()->constrainTranslation(deltaT, this);

  setTranslation(this->translation() + deltaT);
  translation = this->translation();
}

/*! Same as setRotation(), but \p rotation is modified so that the potential
  constraint() of the Frame is satisfied.

  Emits the modified() signal. See also setTranslationWithConstraint() and
  setOrientationWithConstraint(). */
CGAL_INLINE_FUNCTION
void Frame::setRotationWithConstraint(Quaternion &rotation) {
  Quaternion deltaQ = this->rotation().inverse() * rotation;
  if (constraint())
    constraint()->constrainRotation(deltaQ, this);

  // Prevent numerical drift
  deltaQ.normalize();

  setRotation(this->rotation() * deltaQ);
  q_.normalize();
  rotation = this->rotation();
}

/*! Same as setTranslationAndRotation(), but \p translation and \p orientation
  are modified to satisfy the constraint(). Emits the modified() signal. */
CGAL_INLINE_FUNCTION
void Frame::setTranslationAndRotationWithConstraint(Vec &translation,
                                                    Quaternion &rotation) {
  Vec deltaT = translation - this->translation();
  Quaternion deltaQ = this->rotation().inverse() * rotation;

  if (constraint()) {
    constraint()->constrainTranslation(deltaT, this);
    constraint()->constrainRotation(deltaQ, this);
  }

  // Prevent numerical drift
  deltaQ.normalize();

  t_ += deltaT;
  q_ *= deltaQ;
  q_.normalize();

  translation = this->translation();
  rotation = this->rotation();

  Q_EMIT modified();
}

/*! Same as setPosition(), but \p position is modified so that the potential
  constraint() of the Frame is satisfied. See also
  setOrientationWithConstraint() and setTranslationWithConstraint(). */
CGAL_INLINE_FUNCTION
void Frame::setPositionWithConstraint(Vec &position) {
  if (referenceFrame())
    position = referenceFrame()->coordinatesOf(position);

  setTranslationWithConstraint(position);
}

/*! Same as setOrientation(), but \p orientation is modified so that the
  potential constraint() of the Frame is satisfied. See also
  setPositionWithConstraint() and setRotationWithConstraint(). */
CGAL_INLINE_FUNCTION
void Frame::setOrientationWithConstraint(Quaternion &orientation) {
  if (referenceFrame())
    orientation = referenceFrame()->orientation().inverse() * orientation;

  setRotationWithConstraint(orientation);
}

/*! Same as setPositionAndOrientation() but \p position and \p orientation are
modified to satisfy the constraint. Emits the modified() signal. */
CGAL_INLINE_FUNCTION
void Frame::setPositionAndOrientationWithConstraint(Vec &position,
                                                    Quaternion &orientation) {
  if (referenceFrame()) {
    position = referenceFrame()->coordinatesOf(position);
    orientation = referenceFrame()->orientation().inverse() * orientation;
  }
  setTranslationAndRotationWithConstraint(position, orientation);
}

///////////////////////////// REFERENCE FRAMES
//////////////////////////////////////////

/*! Sets the referenceFrame() of the Frame.

The Frame translation() and rotation() are then defined in the referenceFrame()
coordinate system. Use position() and orientation() to express these in the
world coordinate system.

Emits the modified() signal if \p refFrame differs from the current
referenceFrame().

Using this method, you can create a hierarchy of Frames. This hierarchy needs to
be a tree, which root is the world coordinate system (i.e. a \c nullptr
referenceFrame()). A warning is printed and no action is performed if setting \p
refFrame as the referenceFrame() would create a loop in the Frame hierarchy (see
settingAsReferenceFrameWillCreateALoop()). */
CGAL_INLINE_FUNCTION
void Frame::setReferenceFrame(const Frame *const refFrame) {
  if (settingAsReferenceFrameWillCreateALoop(refFrame))
    qWarning("Frame::setReferenceFrame would create a loop in Frame hierarchy");
  else {
    bool identical = (referenceFrame_ == refFrame);
    referenceFrame_ = refFrame;
    if (!identical)
      Q_EMIT modified();
  }
}

/*! Returns \c true if setting \p frame as the Frame's referenceFrame() would
  create a loop in the Frame hierarchy. */
CGAL_INLINE_FUNCTION
bool Frame::settingAsReferenceFrameWillCreateALoop(const Frame *const frame) {
  const Frame *f = frame;
  while (f != nullptr) {
    if (f == this)
      return true;
    f = f->referenceFrame();
  }
  return false;
}

///////////////////////// FRAME TRANSFORMATIONS OF 3D POINTS
/////////////////////////////////

/*! Returns the Frame coordinates of a point \p src defined in the world
 coordinate system (converts from world to Frame).

 inverseCoordinatesOf() performs the inverse conversion. transformOf() converts
 3D vectors instead of 3D coordinates.

 See the <a href="../examples/frameTransform.html">frameTransform example</a>
 for an illustration. */
CGAL_INLINE_FUNCTION
Vec Frame::coordinatesOf(const Vec &src) const {
  if (referenceFrame())
    return localCoordinatesOf(referenceFrame()->coordinatesOf(src));
  else
    return localCoordinatesOf(src);
}

/*! Returns the world coordinates of the point whose position in the Frame
  coordinate system is \p src (converts from Frame to world).

  coordinatesOf() performs the inverse conversion. Use inverseTransformOf() to
  transform 3D vectors instead of 3D coordinates. */
CGAL_INLINE_FUNCTION
Vec Frame::inverseCoordinatesOf(const Vec &src) const {
  const Frame *fr = this;
  Vec res = src;
  while (fr != nullptr) {
    res = fr->localInverseCoordinatesOf(res);
    fr = fr->referenceFrame();
  }
  return res;
}

/*! Returns the Frame coordinates of a point \p src defined in the
  referenceFrame() coordinate system (converts from referenceFrame() to Frame).

  localInverseCoordinatesOf() performs the inverse conversion. See also
  localTransformOf(). */
CGAL_INLINE_FUNCTION
Vec Frame::localCoordinatesOf(const Vec &src) const {
  return rotation().inverseRotate(src - translation());
}

/*! Returns the referenceFrame() coordinates of a point \p src defined in the
 Frame coordinate system (converts from Frame to referenceFrame()).

 localCoordinatesOf() performs the inverse conversion. See also
 localInverseTransformOf(). */
CGAL_INLINE_FUNCTION
Vec Frame::localInverseCoordinatesOf(const Vec &src) const {
  return rotation().rotate(src) + translation();
}

/*! Returns the Frame coordinates of the point whose position in the \p from
  coordinate system is \p src (converts from \p from to Frame).

  coordinatesOfIn() performs the inverse transformation. */
CGAL_INLINE_FUNCTION
Vec Frame::coordinatesOfFrom(const Vec &src, const Frame *const from) const {
  if (this == from)
    return src;
  else if (referenceFrame())
    return localCoordinatesOf(referenceFrame()->coordinatesOfFrom(src, from));
  else
    return localCoordinatesOf(from->inverseCoordinatesOf(src));
}

/*! Returns the \p in coordinates of the point whose position in the Frame
  coordinate system is \p src (converts from Frame to \p in).

  coordinatesOfFrom() performs the inverse transformation. */
CGAL_INLINE_FUNCTION
Vec Frame::coordinatesOfIn(const Vec &src, const Frame *const in) const {
  const Frame *fr = this;
  Vec res = src;
  while ((fr != nullptr) && (fr != in)) {
    res = fr->localInverseCoordinatesOf(res);
    fr = fr->referenceFrame();
  }

  if (fr != in)
    // in was not found in the branch of this, res is now expressed in the world
    // coordinate system. Simply convert to in coordinate system.
    res = in->coordinatesOf(res);

  return res;
}

////// qreal[3] versions

/*! Same as coordinatesOf(), but with \c qreal parameters. */
CGAL_INLINE_FUNCTION
void Frame::getCoordinatesOf(const qreal src[3], qreal res[3]) const {
  const Vec r = coordinatesOf(Vec(src));
  for (int i = 0; i < 3; ++i)
    res[i] = r[i];
}

/*! Same as inverseCoordinatesOf(), but with \c qreal parameters. */
CGAL_INLINE_FUNCTION
void Frame::getInverseCoordinatesOf(const qreal src[3], qreal res[3]) const {
  const Vec r = inverseCoordinatesOf(Vec(src));
  for (int i = 0; i < 3; ++i)
    res[i] = r[i];
}

/*! Same as localCoordinatesOf(), but with \c qreal parameters. */
CGAL_INLINE_FUNCTION
void Frame::getLocalCoordinatesOf(const qreal src[3], qreal res[3]) const {
  const Vec r = localCoordinatesOf(Vec(src));
  for (int i = 0; i < 3; ++i)
    res[i] = r[i];
}

/*! Same as localInverseCoordinatesOf(), but with \c qreal parameters. */
CGAL_INLINE_FUNCTION
void Frame::getLocalInverseCoordinatesOf(const qreal src[3],
                                         qreal res[3]) const {
  const Vec r = localInverseCoordinatesOf(Vec(src));
  for (int i = 0; i < 3; ++i)
    res[i] = r[i];
}

/*! Same as coordinatesOfIn(), but with \c qreal parameters. */
CGAL_INLINE_FUNCTION
void Frame::getCoordinatesOfIn(const qreal src[3], qreal res[3],
                               const Frame *const in) const {
  const Vec r = coordinatesOfIn(Vec(src), in);
  for (int i = 0; i < 3; ++i)
    res[i] = r[i];
}

/*! Same as coordinatesOfFrom(), but with \c qreal parameters. */
CGAL_INLINE_FUNCTION
void Frame::getCoordinatesOfFrom(const qreal src[3], qreal res[3],
                                 const Frame *const from) const {
  const Vec r = coordinatesOfFrom(Vec(src), from);
  for (int i = 0; i < 3; ++i)
    res[i] = r[i];
}

///////////////////////// FRAME TRANSFORMATIONS OF VECTORS
/////////////////////////////////

/*! Returns the Frame transform of a vector \p src defined in the world
 coordinate system (converts vectors from world to Frame).

 inverseTransformOf() performs the inverse transformation. coordinatesOf()
 converts 3D coordinates instead of 3D vectors (here only the rotational part of
 the transformation is taken into account).

 See the <a href="../examples/frameTransform.html">frameTransform example</a>
 for an illustration. */
CGAL_INLINE_FUNCTION
Vec Frame::transformOf(const Vec &src) const {
  if (referenceFrame())
    return localTransformOf(referenceFrame()->transformOf(src));
  else
    return localTransformOf(src);
}

/*! Returns the world transform of the vector whose coordinates in the Frame
  coordinate system is \p src (converts vectors from Frame to world).

  transformOf() performs the inverse transformation. Use inverseCoordinatesOf()
  to transform 3D coordinates instead of 3D vectors. */
CGAL_INLINE_FUNCTION
Vec Frame::inverseTransformOf(const Vec &src) const {
  const Frame *fr = this;
  Vec res = src;
  while (fr != nullptr) {
    res = fr->localInverseTransformOf(res);
    fr = fr->referenceFrame();
  }
  return res;
}

/*! Returns the Frame transform of a vector \p src defined in the
  referenceFrame() coordinate system (converts vectors from referenceFrame() to
  Frame).

  localInverseTransformOf() performs the inverse transformation. See also
  localCoordinatesOf(). */
CGAL_INLINE_FUNCTION
Vec Frame::localTransformOf(const Vec &src) const {
  return rotation().inverseRotate(src);
}

/*! Returns the referenceFrame() transform of a vector \p src defined in the
 Frame coordinate system (converts vectors from Frame to referenceFrame()).

 localTransformOf() performs the inverse transformation. See also
 localInverseCoordinatesOf(). */
CGAL_INLINE_FUNCTION
Vec Frame::localInverseTransformOf(const Vec &src) const {
  return rotation().rotate(src);
}

/*! Returns the Frame transform of the vector whose coordinates in the \p from
  coordinate system is \p src (converts vectors from \p from to Frame).

  transformOfIn() performs the inverse transformation. */
CGAL_INLINE_FUNCTION
Vec Frame::transformOfFrom(const Vec &src, const Frame *const from) const {
  if (this == from)
    return src;
  else if (referenceFrame())
    return localTransformOf(referenceFrame()->transformOfFrom(src, from));
  else
    return localTransformOf(from->inverseTransformOf(src));
}

/*! Returns the \p in transform of the vector whose coordinates in the Frame
  coordinate system is \p src (converts vectors from Frame to \p in).

  transformOfFrom() performs the inverse transformation. */
CGAL_INLINE_FUNCTION
Vec Frame::transformOfIn(const Vec &src, const Frame *const in) const {
  const Frame *fr = this;
  Vec res = src;
  while ((fr != nullptr) && (fr != in)) {
    res = fr->localInverseTransformOf(res);
    fr = fr->referenceFrame();
  }

  if (fr != in)
    // in was not found in the branch of this, res is now expressed in the world
    // coordinate system. Simply convert to in coordinate system.
    res = in->transformOf(res);

  return res;
}

/////////////////  qreal[3] versions  //////////////////////

/*! Same as transformOf(), but with \c qreal parameters. */
CGAL_INLINE_FUNCTION
void Frame::getTransformOf(const qreal src[3], qreal res[3]) const {
  Vec r = transformOf(Vec(src));
  for (int i = 0; i < 3; ++i)
    res[i] = r[i];
}

/*! Same as inverseTransformOf(), but with \c qreal parameters. */
CGAL_INLINE_FUNCTION
void Frame::getInverseTransformOf(const qreal src[3], qreal res[3]) const {
  Vec r = inverseTransformOf(Vec(src));
  for (int i = 0; i < 3; ++i)
    res[i] = r[i];
}

/*! Same as localTransformOf(), but with \c qreal parameters. */
CGAL_INLINE_FUNCTION
void Frame::getLocalTransformOf(const qreal src[3], qreal res[3]) const {
  Vec r = localTransformOf(Vec(src));
  for (int i = 0; i < 3; ++i)
    res[i] = r[i];
}

/*! Same as localInverseTransformOf(), but with \c qreal parameters. */
CGAL_INLINE_FUNCTION
void Frame::getLocalInverseTransformOf(const qreal src[3], qreal res[3]) const {
  Vec r = localInverseTransformOf(Vec(src));
  for (int i = 0; i < 3; ++i)
    res[i] = r[i];
}

/*! Same as transformOfIn(), but with \c qreal parameters. */
CGAL_INLINE_FUNCTION
void Frame::getTransformOfIn(const qreal src[3], qreal res[3],
                             const Frame *const in) const {
  Vec r = transformOfIn(Vec(src), in);
  for (int i = 0; i < 3; ++i)
    res[i] = r[i];
}

/*! Same as transformOfFrom(), but with \c qreal parameters. */
CGAL_INLINE_FUNCTION
void Frame::getTransformOfFrom(const qreal src[3], qreal res[3],
                               const Frame *const from) const {
  Vec r = transformOfFrom(Vec(src), from);
  for (int i = 0; i < 3; ++i)
    res[i] = r[i];
}


/////////////////////////////////   ALIGN   /////////////////////////////////

/*! Aligns the Frame with \p frame, so that two of their axis are parallel.

If one of the X, Y and Z axis of the Frame is almost parallel to any of the X,
Y, or Z axis of \p frame, the Frame is rotated so that these two axis actually
become parallel.

If, after this first rotation, two other axis are also almost parallel, a second
alignment is performed. The two frames then have identical orientations, up to
90 degrees rotations.

\p threshold measures how close two axis must be to be considered parallel. It
is compared with the absolute values of the dot product of the normalized axis.
As a result, useful range is sqrt(2)/2 (systematic alignment) to 1 (no
alignment).

When \p move is set to \c true, the Frame position() is also affected by the
alignment. The new Frame's position() is such that the \p frame position
(computed with coordinatesOf(), in the Frame coordinates system) does not
change.

\p frame may be \c nullptr and then represents the world coordinate system (same
convention than for the referenceFrame()).

The rotation (and translation when \p move is \c true) applied to the Frame are
filtered by the possible constraint(). */
CGAL_INLINE_FUNCTION
void Frame::alignWithFrame(const Frame *const frame, bool move,
                           qreal threshold) {
  Vec directions[2][3];
  for (unsigned short d = 0; d < 3; ++d) {
    Vec dir((d == 0) ? 1.0 : 0.0, (d == 1) ? 1.0 : 0.0, (d == 2) ? 1.0 : 0.0);
    if (frame)
      directions[0][d] = frame->inverseTransformOf(dir);
    else
      directions[0][d] = dir;
    directions[1][d] = inverseTransformOf(dir);
  }

  qreal maxProj = 0.0;
  qreal proj;
  unsigned short index[2];
  index[0] = index[1] = 0;
  for (unsigned short i = 0; i < 3; ++i)
    for (unsigned short j = 0; j < 3; ++j)
      if ((proj = fabs(directions[0][i] * directions[1][j])) >= maxProj) {
        index[0] = i;
        index[1] = j;
        maxProj = proj;
      }

  Frame old;
  old = *this;

  qreal coef = directions[0][index[0]] * directions[1][index[1]];
  if (fabs(coef) >= threshold) {
    const Vec axis = cross(directions[0][index[0]], directions[1][index[1]]);
    qreal angle = asin(axis.norm());
    if (coef >= 0.0)
      angle = -angle;
    rotate(rotation().inverse() * Quaternion(axis, angle) * orientation());

    // Try to align an other axis direction
    unsigned short d = (unsigned short)((index[1] + 1) % 3);
    Vec dir((d == 0) ? 1.0 : 0.0, (d == 1) ? 1.0 : 0.0, (d == 2) ? 1.0 : 0.0);
    dir = inverseTransformOf(dir);

    qreal max = 0.0;
    for (unsigned short i = 0; i < 3; ++i) {
      qreal proj = fabs(directions[0][i] * dir);
      if (proj > max) {
        index[0] = i;
        max = proj;
      }
    }

    if (max >= threshold) {
      const Vec axis = cross(directions[0][index[0]], dir);
      qreal angle = asin(axis.norm());
      if (directions[0][index[0]] * dir >= 0.0)
        angle = -angle;
      rotate(rotation().inverse() * Quaternion(axis, angle) * orientation());
    }
  }

  if (move) {
    Vec center;
    if (frame)
      center = frame->position();

    translate(center - orientation().rotate(old.coordinatesOf(center)) -
              translation());
  }
}

/*! Translates the Frame so that its position() lies on the line defined by \p
origin and \p direction (defined in the world coordinate system).

Simply uses an orthogonal projection. \p direction does not need to be
normalized. */
CGAL_INLINE_FUNCTION
void Frame::projectOnLine(const Vec &origin, const Vec &direction) {
  // If you are trying to find a bug here, because of memory problems, you waste
  // your time. This is a bug in the gcc 3.3 compiler. Compile the library in
  // debug mode and test. Uncommenting this line also seems to solve the
  // problem. Horrible. cout << "position = " << position() << endl; If you
  // found a problem or are using a different compiler, please let me know.
  const Vec shift = origin - position();
  Vec proj = shift;
  proj.projectOnAxis(direction);
  translate(shift - proj);
}
}}
