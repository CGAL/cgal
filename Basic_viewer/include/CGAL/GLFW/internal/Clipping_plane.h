#ifndef CGAL_CLIPPING_PLANE_H
#define CGAL_CLIPPING_PLANE_H

#include "utils.h"

class Clipping_plane 
{
public:
  
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
    m_smoothRotationFactor(0.2f),
    m_smoothTranslationFactor(0.2f),
    m_size(1.0f)
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
  const mat4f& matrix() const
  {
    mat4f translation = transform::translation(m_translation);
    mat3f rotation3x3 = m_orientation.toRotationMatrix();
    mat4f rotation = mat4f::Identity();
    rotation.block<3, 3>(0, 0) = rotation3x3;

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
    m_targetPitch += y;
    m_targetYaw += x;
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
    m_targetTranslation -= m_size * direction * s * m_translationSpeed; 
  }

  inline 
  void translation(const float s) 
  {
    vec3f forward = m_orientation * vec3f::UnitZ();
    m_targetTranslation -= m_size * forward * s * m_translationSpeed; 
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
};

#endif // CGAL_CLIPPING_PLANE_H
