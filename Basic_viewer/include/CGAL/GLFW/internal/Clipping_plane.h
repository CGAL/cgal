#ifndef CGAL_CLIPPING_PLANE_H
#define CGAL_CLIPPING_PLANE_H

#include "utils.h"

class Clipping_plane 
{
public:
  enum class Constraint_axis  
  { 
    NO_CONSTRAINT,
    RIGHT_AXIS,
    UP_AXIS,
  };
public:
  Clipping_plane() : 
    m_orientation(quatf::Identity()), 
    m_translation(vec3f::Zero()),  
    m_targetTranslation(vec3f::Zero()),
    m_upAxis(vec3f::UnitY()),
    m_rightAxis(vec3f::UnitX()),
    m_pitch(0.0f),
    m_targetPitch(0.0f),
    m_yaw(0.0f),
    m_targetYaw(0.0f),
    m_translationSpeed(1.0f),
    m_rotationSpeed(200.0f),
    m_smoothRotationFactor(0.5f),
    m_smoothTranslationFactor(0.5f),
    m_size(1.0f),
    m_constraintAxis(Constraint_axis::NO_CONSTRAINT)
  {
  }

  inline 
  void reset_translation() 
  {
    m_translation = vec3f::Zero();
    m_targetTranslation = vec3f::Zero();
    m_translationSpeed = 1.0f;
  }

  inline 
  void reset_rotation() 
  {
    m_orientation = quatf::Identity(); 
    m_pitch = 0.0f;
    m_targetPitch = 0.0f;
    m_yaw = 0.0f;
    m_targetYaw = 0.0f;
    m_rotationSpeed = 200.0f;

    m_upAxis = vec3f::UnitY();
    m_rightAxis = vec3f::UnitX();
  }

  inline 
  void reset_all()
  {
    reset_translation();
    reset_rotation();
    m_smoothRotationFactor = 0.1f;
    m_smoothTranslationFactor = 0.1f;
    m_size = 1.0f;
  }

  inline 
  mat4f matrix() const
  {
    mat4f translation = transform::translation(m_translation);
    mat4f rotation = transform::rotation(m_orientation);

    return translation * rotation;
  }

  inline 
  void update(const float dt)
  {
    float smoothPitch = m_pitch + m_smoothRotationFactor * (m_targetPitch - m_pitch);   
    float smoothYaw = m_yaw + m_smoothRotationFactor * (m_targetYaw - m_yaw);   

    float pitchDelta = radians((smoothPitch - m_pitch) * m_rotationSpeed) * dt;
    float yawDelta = radians((smoothYaw - m_yaw) * m_rotationSpeed) * dt;

    quatf pitchQuaternion(Eigen::AngleAxisf(pitchDelta, m_rightAxis));
    quatf yawQuaternion(Eigen::AngleAxisf(yawDelta, m_upAxis));
    m_orientation = pitchQuaternion * yawQuaternion * m_orientation;

    m_pitch = smoothPitch;
    m_yaw = smoothYaw;

    m_translation += m_smoothTranslationFactor * (m_targetTranslation - m_translation);
  }

  inline 
  void rotation(const float x, const float y) 
  {
    if (
      m_constraintAxis == Constraint_axis::NO_CONSTRAINT ||
      m_constraintAxis == Constraint_axis::RIGHT_AXIS) 
    {
      m_targetPitch += y;
    }
    
    if (
      m_constraintAxis == Constraint_axis::NO_CONSTRAINT ||
      m_constraintAxis == Constraint_axis::UP_AXIS)
    {
      m_targetYaw += x;
    }
  }

  inline 
  void translation(const float x, const float y) 
  {
    m_targetTranslation -= m_size * m_rightAxis * x * m_translationSpeed; 
    m_targetTranslation -= m_size * m_upAxis * y * m_translationSpeed; 
  }

  inline 
  void translation(const vec3f& direction, const float s) 
  {
    m_targetTranslation -= m_size * direction.normalized() * s * m_translationSpeed; 
  }

  inline 
  void translation(const float s) 
  {
    vec3f forward = (m_orientation * vec3f::UnitZ()).normalized();
    m_targetTranslation -= m_size * forward * s * m_translationSpeed; 
  }

  inline 
  void switch_constraint_axis() 
  {
    switch(m_constraintAxis)
    {
      case Constraint_axis::NO_CONSTRAINT: 
        m_constraintAxis = Constraint_axis::RIGHT_AXIS;
        break;
      case Constraint_axis::RIGHT_AXIS: 
        m_constraintAxis = Constraint_axis::UP_AXIS;
        break;
      case Constraint_axis::UP_AXIS: 
        m_constraintAxis = Constraint_axis::NO_CONSTRAINT;
        break;
    }
  }

  inline void set_size(const float size) { m_size = size; }
  inline void set_up_axis(const vec3f& upAxis) { m_upAxis = upAxis; }
  inline void set_right_axis(const vec3f& rightAxis) { m_rightAxis = rightAxis; }

private:
  quatf m_orientation;

  vec3f m_translation;
  vec3f m_targetTranslation;

  vec3f m_upAxis;
  vec3f m_rightAxis;

  float m_pitch;
  float m_targetPitch;
  float m_yaw;
  float m_targetYaw;

  float m_rotationSpeed;
  float m_translationSpeed;

  float m_smoothTranslationFactor;
  float m_smoothRotationFactor;

  float m_size; 

  Constraint_axis m_constraintAxis;
};

#endif // CGAL_CLIPPING_PLANE_H
