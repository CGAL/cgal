/****************************************************************************

 Copyright (c) 2018  GeometryFactory Sarl (France).
 Copyright (C) 2002-2014 Gilles Debunne. All rights reserved.

 This file is part of a fork of the QGLViewer library version 2.7.0.

*****************************************************************************/
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-only
#ifndef QGLVIEWER_KEY_FRAME_INTERPOLATOR_H
#define QGLVIEWER_KEY_FRAME_INTERPOLATOR_H

#include <QObject>
#include <QTimer>
#include <CGAL/export/Qt.h>
#include <CGAL/Qt/quaternion.h>
// Not actually needed, but some bad compilers (Microsoft VS6) complain.
//#include <CGAL/Qt/frame.h>

// If you compiler complains about incomplete type, uncomment the next line
// #include "frame.h"
// and comment "class Frame;" 3 lines below

namespace CGAL{
namespace qglviewer {
class Frame;
/*! \brief A keyFrame Catmull-Rom Frame interpolator.
  \class KeyFrameInterpolator keyFrameInterpolator.h
  CGAL::QGLViewer/keyFrameInterpolator.h

  A KeyFrameInterpolator holds keyFrames (that define a path) and a pointer to a
  Frame of your application (which will be interpolated). When the user
  startInterpolation(), the KeyFrameInterpolator regularly updates the frame()
  position and orientation along the path.

  Here is a typical utilization example (see also the <a
  href="../examples/keyFrames.html">keyFrames example</a>): \code


  init()
  {
        // The KeyFrameInterpolator kfi is given the Frame that it will drive
  over time. kfi = new KeyFrameInterpolator( new Frame() ); kfi->addKeyFrame(
  Frame( Vec(1,0,0), Quaternion() ) ); kfi->addKeyFrame( new Frame( Vec(2,1,0),
  Quaternion() ) );
        // ...and so on for all the keyFrames.

        // Ask for a display update after each update of the
  KeyFrameInterpolator connect(kfi, SIGNAL(interpolated()), SLOT(update()));

        kfi->startInterpolation();
  }

  draw()
  {
        glPushMatrix();
        glMultMatrixd( kfi->frame()->matrix() );
        // Draw your object here. Its position and orientation are interpolated.
        glPopMatrix();
  }
  \endcode

  The keyFrames are defined by a Frame and a time, expressed in seconds. The
  Frame can be provided as a const reference or as a pointer to a Frame (see the
  addKeyFrame() methods). In the latter case, the path will automatically be
  updated when the Frame is modified (using the Frame::modified() signal).

  The time has to be monotonously increasing over keyFrames. When
  interpolationSpeed() equals 1.0 (default value), these times correspond to
  actual user's seconds during interpolation (provided that your main loop is
  fast enough). The interpolation is then real-time: the keyFrames will be
  reached at their keyFrameTime().

  <h3>Interpolation details</h3>

  When the user startInterpolation(), a timer is started which will update the
  frame()'s position and orientation every interpolationPeriod() milliseconds.
  This update increases the interpolationTime() by interpolationPeriod() *
  interpolationSpeed() milliseconds.

  Note that this mechanism ensures that the number of interpolation steps is
  constant and equal to the total path duration() divided by the
  interpolationPeriod() * interpolationSpeed(). This is especially useful for
  benchmarking or movie creation (constant number of snapshots).

  During the interpolation, the KeyFrameInterpolator emits an interpolated()
  signal, which will usually be connected to the CGAL::QGLViewer::update() slot. The
  interpolation is stopped when interpolationTime() is greater than the
  lastTime() (unless loopInterpolation() is \c true) and the endReached() signal
  is then emitted.

  Note that a Camera has Camera::keyFrameInterpolator(), that can be used to
  drive the Camera along a path, or to restore a saved position (a path made of
  a single keyFrame). Press Alt+Fx to define a new keyFrame for path x. Pressing
  Fx plays/pauses path interpolation. See CGAL::QGLViewer::pathKey() and the <a
  href="../keyboard.html">keyboard page</a> for details.

  \attention If a Constraint is attached to the frame() (see
  Frame::constraint()), it should be deactivated before
  interpolationIsStarted(), otherwise the interpolated motion (computed as if
  there was no constraint) will probably be erroneous.

  <h3>Retrieving interpolated values</h3>

  This code defines a KeyFrameInterpolator, and displays the positions that will
  be followed by the frame() along the path: \code KeyFrameInterpolator kfi( new
  Frame() );
  // calls to kfi.addKeyFrame() to define the path.

  const qreal deltaTime = 0.04; // output a position every deltaTime seconds
  for (qreal time=kfi.firstTime(); time<=kfi.lastTime(); time += deltaTime)
  {
        kfi.interpolateAtTime(time);
        cout << "t=" << time << "\tpos=" << kfi.frame()->position() << endl;
  }
  \endcode
  You may want to temporally disconnect the \c kfi interpolated() signal from
  the CGAL::QGLViewer::update() slot before calling this code. \nosubgrouping */
class CGAL_QT_EXPORT KeyFrameInterpolator : public QObject {
  // todo closedPath, insertKeyFrames, deleteKeyFrame, replaceKeyFrame
  Q_OBJECT

public:
  KeyFrameInterpolator(Frame *fr = nullptr);
  virtual ~KeyFrameInterpolator();

Q_SIGNALS:
  /*! This signal is emitted whenever the frame() state is interpolated.

  The emission of this signal triggers the synchronous emission of the frame()
  Frame::interpolated() signal, which may also be useful.

  This signal should especially be connected to your CGAL::QGLViewer::update() slot,
  so that the display is updated after every update of the KeyFrameInterpolator
  frame(): \code connect(myKeyFrameInterpolator, SIGNAL(interpolated()),
  SLOT(update())); \endcode Use the CGAL::QGLViewer::QGLViewerPool() to connect the
  signal to all the viewers.

  Note that the CGAL::QGLViewer::camera() Camera::keyFrameInterpolator() created using
  CGAL::QGLViewer::pathKey() have their interpolated() signals automatically connected
  to the CGAL::QGLViewer::update() slot. */
  void interpolated();

  /*! This signal is emitted when the interpolation reaches the first (when
  interpolationSpeed() is negative) or the last keyFrame.

  When loopInterpolation() is \c true, interpolationTime() is reset and the
  interpolation continues. It otherwise stops. */
  void endReached();

  /*! @name Path creation */
  //@{
public Q_SLOTS:
  void addKeyFrame(const Frame &frame);
  void addKeyFrame(const Frame &frame, qreal time);

  void addKeyFrame(const Frame *const frame);
  void addKeyFrame(const Frame *const frame, qreal time);

  void deletePath();
  //@}

  /*! @name Associated Frame */
  //@{
public:
  /*! Returns the associated Frame and that is interpolated by the
  KeyFrameInterpolator.

  When interpolationIsStarted(), this Frame's position and orientation will
  regularly be updated by a timer, so that they follow the KeyFrameInterpolator
  path.

  Set using setFrame() or with the KeyFrameInterpolator constructor. */
  Frame *frame() const { return frame_; }

public Q_SLOTS:
  void setFrame(Frame *const frame);
  //@}

  /*! @name Path parameters */
  //@{
public:
  Frame keyFrame(int index) const;
  qreal keyFrameTime(int index) const;
  /*! Returns the number of keyFrames used by the interpolation. Use
   * addKeyFrame() to add new keyFrames. */
  int numberOfKeyFrames() const { return keyFrame_.count(); }
  qreal duration() const;
  qreal firstTime() const;
  qreal lastTime() const;
  //@}

  /*! @name Interpolation parameters */
  //@{
public:
  /*! Returns the current interpolation time (in seconds) along the
  KeyFrameInterpolator path.

  This time is regularly updated when interpolationIsStarted(). Can be set
  directly with setInterpolationTime() or interpolateAtTime(). */
  qreal interpolationTime() const { return interpolationTime_; }
  /*! Returns the current interpolation speed.

  Default value is 1.0, which means keyFrameTime() will be matched during the
  interpolation (provided that your main loop is fast enough).

  A negative value will result in a reverse interpolation of the keyFrames. See
  also interpolationPeriod(). */
  qreal interpolationSpeed() const { return interpolationSpeed_; }
  /*! Returns the current interpolation period, expressed in milliseconds.

  The update of the frame() state will be done by a timer at this period when
  interpolationIsStarted().

  This period (multiplied by interpolationSpeed()) is added to the
  interpolationTime() at each update, and the frame() state is modified
  accordingly (see interpolateAtTime()). Default value is 40 milliseconds. */
  int interpolationPeriod() const { return period_; }
  /*! Returns \c true when the interpolation is played in an infinite loop.

  When \c false (default), the interpolation stops when interpolationTime()
  reaches firstTime() (with negative interpolationSpeed()) or lastTime().

  interpolationTime() is otherwise reset to firstTime() (+ interpolationTime() -
  lastTime()) (and inversely for negative interpolationSpeed()) and
  interpolation continues.

  In both cases, the endReached() signal is emitted. */
  bool loopInterpolation() const { return loopInterpolation_; }
#ifndef DOXYGEN
  /*! Whether or not (default) the path defined by the keyFrames is a closed
  loop. When \c true, the last and the first KeyFrame are linked by a new spline
  segment.

  Use setLoopInterpolation() to create a continuous animation over the entire
  path. \attention The closed path feature is not yet implemented. */
  bool closedPath() const { return closedPath_; }
#endif
public Q_SLOTS:
  /*! Sets the interpolationTime().

  \attention The frame() state is not affected by this method. Use this function
  to define the starting time of a future interpolation (see
  startInterpolation()). Use interpolateAtTime() to actually interpolate at a
  given time. */
  void setInterpolationTime(qreal time) { interpolationTime_ = time; }
  /*! Sets the interpolationSpeed(). Negative or null values are allowed. */
  void setInterpolationSpeed(qreal speed) { interpolationSpeed_ = speed; }
  /*! Sets the interpolationPeriod(). */
  void setInterpolationPeriod(int period) { period_ = period; }
  /*! Sets the loopInterpolation() value. */
  void setLoopInterpolation(bool loop = true) { loopInterpolation_ = loop; }
#ifndef DOXYGEN
  /*! Sets the closedPath() value. \attention The closed path feature is not yet
   * implemented. */
  void setClosedPath(bool closed = true) { closedPath_ = closed; }
#endif
  //@}

  /*! @name Interpolation */
  //@{
public:
  /*! Returns \c true when the interpolation is being performed. Use
  startInterpolation(), stopInterpolation() or toggleInterpolation() to modify
  this state. */
  bool interpolationIsStarted() const { return interpolationStarted_; }
public Q_SLOTS:
  void startInterpolation(int period = -1);
  void stopInterpolation();
  void resetInterpolation();
  /*! Calls startInterpolation() or stopInterpolation(), depending on
   * interpolationIsStarted(). */
  void toggleInterpolation() {
    if (interpolationIsStarted())
      stopInterpolation();
    else
      startInterpolation();
  }
  virtual void interpolateAtTime(qreal time);
  //@}

private Q_SLOTS:
  virtual void update();
  virtual void invalidateValues() {
    valuesAreValid_ = false;
    pathIsValid_ = false;
    splineCacheIsValid_ = false;
  }

private:
  // Copy constructor and opertor= are declared private and undefined
  // Prevents everyone from trying to use them
  // KeyFrameInterpolator(const KeyFrameInterpolator& kfi);
  // KeyFrameInterpolator& operator=(const KeyFrameInterpolator& kfi);

  void updateCurrentKeyFrameForTime(qreal time);
  void updateModifiedFrameValues();
  void updateSplineCache();

#ifndef DOXYGEN
  // Internal private KeyFrame representation
  class KeyFrame {
  public:
    KeyFrame(const Frame &fr, qreal t);
    KeyFrame(const Frame *fr, qreal t);

    Vec position() const { return p_; }
    Quaternion orientation() const { return q_; }
    Vec tgP() const { return tgP_; }
    Quaternion tgQ() const { return tgQ_; }
    qreal time() const { return time_; }
    const Frame *frame() const { return frame_; }
    void updateValuesFromPointer();
    void flipOrientationIfNeeded(const Quaternion &prev);
    void computeTangent(const KeyFrame *const prev, const KeyFrame *const next);

  private:
    Vec p_, tgP_;
    Quaternion q_, tgQ_;
    qreal time_;
    const Frame *const frame_;
  };
#endif

  // K e y F r a m e s
  mutable QList<KeyFrame *> keyFrame_;
  QMutableListIterator<KeyFrame *> *currentFrame_[4];
  QList<Frame> path_;

  // A s s o c i a t e d   f r a m e
  Frame *frame_;

  // R h y t h m
  QTimer timer_;
  int period_;
  qreal interpolationTime_;
  qreal interpolationSpeed_;
  bool interpolationStarted_;

  // M i s c
  bool closedPath_;
  bool loopInterpolation_;

  // C a c h e d   v a l u e s   a n d   f l a g s
  bool pathIsValid_;
  bool valuesAreValid_;
  bool currentFrameValid_;
  bool splineCacheIsValid_;
  Vec v1, v2;
};

}} // namespace CGAL::qglviewer

#endif // QGLVIEWER_KEY_FRAME_INTERPOLATOR_H
