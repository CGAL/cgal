/****************************************************************************

 Copyright (c) 2018  GeometryFactory Sarl (France).
 Copyright (C) 2002-2014 Gilles Debunne. All rights reserved.

 This file is part of a fork of the QGLViewer library version 2.7.0.

*****************************************************************************/
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-only
#ifndef QGLVIEWER_CONSTRAINT_H
#define QGLVIEWER_CONSTRAINT_H

#include <CGAL/Qt/quaternion.h>
#include <CGAL/Qt/vec.h>
#include <CGAL/export/Qt.h>

namespace CGAL{
namespace qglviewer {
class Frame;
class Camera;
/*! \brief An interface class for Frame constraints.
  \class Constraint constraint.h CGAL::QGLViewer/constraint.h

  This class defines the interface for the Constraints that can be applied to a
  Frame to limit its motion. Use Frame::setConstraint() to associate a
  Constraint to a Frame (default is a \c nullptr Frame::constraint()).

  <h3>How does it work ?</h3>

  The Constraint acts as a filter on the translation and rotation Frame
  increments. constrainTranslation() and constrainRotation() should be
  overloaded to specify the constraint behavior: the desired displacement is
  given as a parameter that can optionally be modified.

  Here is how the Frame::translate() and Frame::rotate() methods use the
  Constraint: \code Frame::translate(Vec& T)
  {
        if (constraint())
          constraint()->constrainTranslation(T, this);
        t += T;
  }

  Frame::rotate(Quaternion& Q)
  {
        if (constraint())
          constraint()->constrainRotation(Q, this);
        q *= Q;
  }
  \endcode

  The default behavior of constrainTranslation() and constrainRotation() is
  empty (meaning no filtering).

  The Frame which uses the Constraint is passed as a parameter to the
  constrainTranslation() and constrainRotation() methods, so that they can have
  access to its current state (mainly Frame::position() and
  Frame::orientation()). It is not \c const for versatility reasons, but
  directly modifying it should be avoided.

  \attention Frame::setTranslation(), Frame::setRotation() and similar methods
  will actually indeed set the frame position and orientation, without taking
  the constraint into account. Use the \e WithConstraint versions of these
  methods to enforce the Constraint.

  <h3>Implemented Constraints</h3>

  Classical axial and plane Constraints are provided for convenience: see the
  LocalConstraint, WorldConstraint and CameraConstraint classes' documentations.

  Try the <a href="../examples/constrainedFrame.html">constrainedFrame</a> and
  <a href="../examples/constrainedCamera.html">constrainedCamera</a> examples
  for an illustration.

  <h3>Creating new Constraints</h3>

  The implementation of a new Constraint class simply consists in overloading
  the filtering methods: \code
  // This Constraint enforces that the Frame cannot have a negative z world
  coordinate. class myConstraint : public Constraint
  {
  public:
        virtual void constrainTranslation(Vec& t, Frame * const fr)
          {
                // Express t in the world coordinate system.
                const Vec tWorld = fr->inverseTransformOf(t);
        if (fr->position().z + tWorld.z < 0.0) // check the new fr z coordinate
          t.z = fr->transformOf(-fr->position().z); // t.z is clamped so that
  next z position is 0.0
          }
  };
  \endcode

  Note that the translation (resp. rotation) parameter passed to
  constrainTranslation() (resp. constrainRotation()) is expressed in the \e
  local Frame coordinate system. Here, we use the Frame::transformOf() and
  Frame::inverseTransformOf() method to convert it to and from the world
  coordinate system.

  Combined constraints can easily be achieved by creating a new class that
  applies the different constraint filters: \code
  myConstraint::constrainTranslation(Vec& v, Frame* const fr)
  {
        constraint1->constrainTranslation(v, fr);
        constraint2->constrainTranslation(v, fr);
        // and so on, with possible branches, tests, loops...
  }
  \endcode
  */
class CGAL_QT_EXPORT Constraint {
public:
  /*! Virtual destructor. Empty. */
  virtual ~Constraint() {}

  /*! Filters the translation applied to the \p frame. This default
  implementation is empty (no filtering).

  Overload this method in your own Constraint class to define a new translation
  constraint. \p frame is the Frame to which is applied the translation. It is
  not defined \c const, but you should refrain from directly changing its value
  in the constraint. Use its Frame::position() and update the \p translation
  accordingly instead.

  \p translation is expressed in local frame coordinate system. Use
  Frame::inverseTransformOf() to express it in the world coordinate system if
  needed. */
  virtual void constrainTranslation(Vec &translation, Frame *const frame) {
    Q_UNUSED(translation);
    Q_UNUSED(frame);
  }
  /*! Filters the rotation applied to the \p frame. This default implementation
  is empty (no filtering).

  Overload this method in your own Constraint class to define a new rotation
  constraint. See constrainTranslation() for details.

  Use Frame::inverseTransformOf() on the \p rotation Quaternion::axis() to
  express \p rotation in the world coordinate system if needed. */
  virtual void constrainRotation(Quaternion &rotation, Frame *const frame) {
    Q_UNUSED(rotation);
    Q_UNUSED(frame);
  }
};

/*!
   \brief An abstract class for Frame Constraints defined by an axis or a plane.
   \class AxisPlaneConstraint constraint.h CGAL::QGLViewer/constraint.h

   AxisPlaneConstraint is an interface for (translation and/or rotation)
   Constraint that are defined by a direction. translationConstraintType() and
   rotationConstraintType() define how this direction should be interpreted: as
   an axis (AxisPlaneConstraint::AXIS) or as a plane normal
   (AxisPlaneConstraint::PLANE). See the Type() documentation for details.

   The three implementations of this class: LocalConstraint, WorldConstraint and
   CameraConstraint differ by the coordinate system in which this direction is
   expressed.

   Different implementations of this class are illustrated in the
   <a href="../examples/constrainedCamera.html">contrainedCamera</a> and
   <a href="../examples/constrainedFrame.html">constrainedFrame</a> examples.

   \attention When applied, the rotational Constraint may not intuitively follow
   the mouse displacement. A solution would be to directly measure the rotation
   angle in screen coordinates, but that would imply to know the
   CGAL::QGLViewer::camera(), so that we can compute the projected coordinates of the
   rotation center (as is done with the CGAL::QGLViewer::SCREEN_ROTATE binding).
   However, adding an extra pointer to the CGAL::QGLViewer::camera() in all the
   AxisPlaneConstraint derived classes (which the user would have to update in a
   multi-viewer application) was judged as an overkill. */
class CGAL_QT_EXPORT AxisPlaneConstraint : public Constraint {
public:
  AxisPlaneConstraint();
  /*! Virtual destructor. Empty. */
  virtual ~AxisPlaneConstraint() {}

  /*! Type lists the different types of translation and rotation constraints
  that are available.

  It specifies the meaning of the constraint direction (see
  translationConstraintDirection() and rotationConstraintDirection()): as an
  axis direction (AxisPlaneConstraint::AXIS) or a plane normal
  (AxisPlaneConstraint::PLANE). AxisPlaneConstraint::FREE means no constraint
  while AxisPlaneConstraint::FORBIDDEN completely forbids the translation and/or
  the rotation.

  See translationConstraintType() and rotationConstraintType().

  \attention The AxisPlaneConstraint::PLANE Type is not valid for rotational
  constraint.

  New derived classes can use their own extended enum for specific constraints:
  \code
  class MyAxisPlaneConstraint : public AxisPlaneConstraint
  {
  public:
    enum MyType { FREE, AXIS, PLANE, FORBIDDEN, CUSTOM };
    virtual void constrainTranslation(Vec &translation, Frame *const frame)
    {
          // translationConstraintType() is simply an int. CUSTOM Type is
  handled seamlessly. switch (translationConstraintType())
          {
            case MyAxisPlaneConstraint::FREE: ... break;
            case MyAxisPlaneConstraint::CUSTOM: ... break;
          }
    };

    MyAxisPlaneConstraint* c = new MyAxisPlaneConstraint();
    // Note the Type conversion
    c->setTranslationConstraintType(AxisPlaneConstraint::Type(MyAxisPlaneConstraint::CUSTOM));
  };
  \endcode */
  enum Type { FREE, AXIS, PLANE, FORBIDDEN };

  /*! @name Translation constraint */
  //@{
  /*! Overloading of Constraint::constrainTranslation(). Empty */
  virtual void constrainTranslation(Vec &translation, Frame *const frame) {
    Q_UNUSED(translation);
    Q_UNUSED(frame);
  };

  void setTranslationConstraint(Type type, const Vec &direction);
  /*! Sets the Type() of the translationConstraintType(). Default is
   * AxisPlaneConstraint::FREE. */
  void setTranslationConstraintType(Type type) {
    translationConstraintType_ = type;
  };
  void setTranslationConstraintDirection(const Vec &direction);

  /*! Returns the translation constraint Type().

  Depending on this value, the Frame will freely translate
  (AxisPlaneConstraint::FREE), will only be able to translate along an axis
  direction (AxisPlaneConstraint::AXIS), will be forced to stay into a plane
  (AxisPlaneConstraint::PLANE) or will not able to translate at all
  (AxisPlaneConstraint::FORBIDDEN).

  Use Frame::setPosition() to define the position of the constrained Frame
  before it gets constrained. */
  Type translationConstraintType() const { return translationConstraintType_; };
  /*! Returns the direction used by the translation constraint.

  It represents the axis direction (AxisPlaneConstraint::AXIS) or the plane
  normal (AxisPlaneConstraint::PLANE) depending on the
  translationConstraintType(). It is undefined for AxisPlaneConstraint::FREE or
  AxisPlaneConstraint::FORBIDDEN.

  The AxisPlaneConstraint derived classes express this direction in different
  coordinate system (camera for CameraConstraint, local for LocalConstraint, and
  world for WorldConstraint). This value can be modified with
  setTranslationConstraintDirection(). */
  Vec translationConstraintDirection() const {
    return translationConstraintDir_;
  };
  //@}

  /*! @name Rotation constraint */
  //@{
  /*! Overloading of Constraint::constrainRotation(). Empty. */
  virtual void constrainRotation(Quaternion &rotation, Frame *const frame) {
    Q_UNUSED(rotation);
    Q_UNUSED(frame);
  };

  void setRotationConstraint(Type type, const Vec &direction);
  void setRotationConstraintType(Type type);
  void setRotationConstraintDirection(const Vec &direction);

  /*! Returns the rotation constraint Type(). */
  Type rotationConstraintType() const { return rotationConstraintType_; };
  /*! Returns the axis direction used by the rotation constraint.

  This direction is defined only when rotationConstraintType() is
  AxisPlaneConstraint::AXIS.

  The AxisPlaneConstraint derived classes express this direction in different
  coordinate system (camera for CameraConstraint, local for LocalConstraint, and
  world for WorldConstraint). This value can be modified with
  setRotationConstraintDirection(). */
  Vec rotationConstraintDirection() const { return rotationConstraintDir_; };
  //@}

private:
  // int and not Type to allow for overloading and new types definition.
  Type translationConstraintType_;
  Type rotationConstraintType_;

  Vec translationConstraintDir_;
  Vec rotationConstraintDir_;
};

/*! \brief An AxisPlaneConstraint defined in the Frame local coordinate system.
  \class LocalConstraint constraint.h CGAL::QGLViewer/constraint.h

  The translationConstraintDirection() and rotationConstraintDirection() are
  expressed in the Frame local coordinate system (see Frame::referenceFrame()).

  See the <a href="../examples/constrainedFrame.html">constrainedFrame</a>
  example for an illustration. */
class CGAL_QT_EXPORT LocalConstraint : public AxisPlaneConstraint {
public:
  /*! Virtual destructor. Empty. */
  virtual ~LocalConstraint(){};

  virtual void constrainTranslation(Vec &translation, Frame *const frame);
  virtual void constrainRotation(Quaternion &rotation, Frame *const frame);
};

/*! \brief An AxisPlaneConstraint defined in the world coordinate system.
        \class WorldConstraint constraint.h CGAL::QGLViewer/constraint.h

  The translationConstraintDirection() and rotationConstraintDirection() are
  expressed in world coordinate system.

  See the <a href="../examples/constrainedFrame.html">constrainedFrame</a> and
  <a href="../examples/multiView.html">multiView</a> examples for an
  illustration. */
class CGAL_QT_EXPORT WorldConstraint : public AxisPlaneConstraint {
public:
  /*! Virtual destructor. Empty. */
  virtual ~WorldConstraint(){};

  virtual void constrainTranslation(Vec &translation, Frame *const frame);
  virtual void constrainRotation(Quaternion &rotation, Frame *const frame);
};

/*! \brief An AxisPlaneConstraint defined in the camera coordinate system.
  \class CameraConstraint constraint.h CGAL::QGLViewer/constraint.h

  The translationConstraintDirection() and rotationConstraintDirection() are
  expressed in the associated camera() coordinate system.

  See the <a href="../examples/constrainedFrame.html">constrainedFrame</a> and
  <a href="../examples/constrainedCamera.html">constrainedCamera</a> examples
  for an illustration. */
class CGAL_QT_EXPORT CameraConstraint : public AxisPlaneConstraint {
public:
  explicit CameraConstraint(const Camera *const camera);
  /*! Virtual destructor. Empty. */
  virtual ~CameraConstraint(){}

  virtual void constrainTranslation(Vec &translation, Frame *const frame);
  virtual void constrainRotation(Quaternion &rotation, Frame *const frame);

  /*! Returns the associated Camera. Set using the CameraConstraint constructor.
   */
  const Camera *camera() const { return camera_; }

private:
  const Camera *const camera_;
};

}} // namespace CGAL::qglviewer

#ifdef CGAL_HEADER_ONLY
//#include <CGAL/Qt/qglviewer_impl_list.h>
#endif // CGAL_HEADER_ONLY
#endif // QGLVIEWER_CONSTRAINT_H
