#ifndef CGAL_CLIPPING_PLANE_H
#define CGAL_CLIPPING_PLANE_H

#include "utils.h"
#include "../bv_settings.h"

class Clipping_plane
{
public:
  enum class ConstraintAxis { NO_CONSTRAINT, RIGHT_AXIS, UP_AXIS };

public:
  void update(const float deltaTime);

  void reset_all();
  void reset_position();
  void reset_orientation();

  mat4f get_matrix() const;

  void rotation(const float x, const float y);

  void translation(const float s);
  void translation(const float x, const float y);
  void translation(const vec3f& direction, const float s);

  void switch_constraint_axis();

  inline float get_transparency() const { return m_transparency; }
  inline float get_rotation_speed() const { return m_rotationSpeed; }
  inline float get_translation_speed() const { return m_translationSpeed; }
  inline std::string get_constraint_axis() const { return m_constraintAxis == ConstraintAxis::NO_CONSTRAINT ? "None" : (m_constraintAxis == ConstraintAxis::RIGHT_AXIS ? "Right" : "Up"); }

  inline void set_size(const float size) { m_size = size; }
  inline void set_up_axis(const vec3f& upAxis) { m_upAxis = upAxis; }
  inline void set_right_axis(const vec3f& rightAxis) { m_rightAxis = rightAxis; }

  inline void increase_rotation_speed(const float deltaTime) { m_rotationSpeed = std::min(m_rotationSpeed+100.f*deltaTime, 360.f); }
  inline void decrease_rotation_speed(const float deltaTime) { m_rotationSpeed = std::max(m_rotationSpeed-100.f*deltaTime, 60.f); }

  inline void increase_translation_speed(const float deltaTime) { m_translationSpeed = std::min(m_translationSpeed+deltaTime, 20.f); }
  inline void decrease_translation_speed(const float deltaTime) { m_translationSpeed = std::max(m_translationSpeed-deltaTime, 0.5f); }

  inline void increase_rotation_smoothness(const float deltaTime) { m_rotationSmoothFactor = std::max(m_rotationSmoothFactor-deltaTime, 0.01f); }
  inline void decrease_rotation_smoothness(const float deltaTime) { m_rotationSmoothFactor = std::min(m_rotationSmoothFactor+deltaTime, 1.f); }

  inline void increase_translation_smoothness(const float deltaTime) { m_translationSmoothFactor = std::max(m_translationSmoothFactor-deltaTime, 0.01f); }
  inline void decrease_translation_smoothness(const float deltaTime) { m_translationSmoothFactor = std::min(m_translationSmoothFactor+deltaTime, 1.f); }

private:
  quatf m_orientation { quatf::Identity() };

  vec3f m_position       { vec3f::Zero() };
  vec3f m_targetPosition { vec3f::Zero() };

  vec3f m_upAxis    { vec3f::UnitY() };
  vec3f m_rightAxis { vec3f::UnitX() };

  float m_pitch       { 0.0f };
  float m_targetPitch { 0.0f };
  float m_yaw         { 0.0f };
  float m_targetYaw   { 0.0f };

  float m_rotationSpeed    { CGAL_CLIPPING_PLANE_ROTATION_SPEED };
  float m_translationSpeed { CGAL_CLIPPING_PLANE_TRANSLATION_SPEED };

  float m_rotationSmoothFactor    { CGAL_CLIPPING_PLANE_ROTATION_SMOOTHNESS };
  float m_translationSmoothFactor { CGAL_CLIPPING_PLANE_TRANSLATION_SMOOTHNESS };

  float m_size { 1.0f }; 

  float m_transparency { CGAL_CLIPPING_PLANE_RENDERING_TRANSPARENCY }; // to what extent the transparent part should be rendered;

  ConstraintAxis m_constraintAxis { ConstraintAxis::NO_CONSTRAINT };
};

/********************METHOD DEFINITIONS********************/

inline 
void Clipping_plane::update(const float deltaTime)
{
  float smoothPitch = m_pitch + m_rotationSmoothFactor * (m_targetPitch - m_pitch);   
  float smoothYaw = m_yaw + m_rotationSmoothFactor * (m_targetYaw - m_yaw);   

  float pitchDelta = utils::radians((smoothPitch - m_pitch) * m_rotationSpeed) * deltaTime;
  float yawDelta = utils::radians((smoothYaw - m_yaw) * m_rotationSpeed) * deltaTime;

  quatf pitchQuaternion(Eigen::AngleAxisf(pitchDelta, m_rightAxis));
  quatf yawQuaternion(Eigen::AngleAxisf(yawDelta, m_upAxis));
  m_orientation = pitchQuaternion * yawQuaternion * m_orientation;

  m_pitch = smoothPitch;
  m_yaw = smoothYaw;

  m_position += m_translationSmoothFactor * (m_targetPosition - m_position);
}

inline 
void Clipping_plane::reset_all()
{
  reset_position();
  reset_orientation();
  
  m_rotationSmoothFactor    = CGAL_CLIPPING_PLANE_ROTATION_SMOOTHNESS;
  m_translationSmoothFactor = CGAL_CLIPPING_PLANE_TRANSLATION_SMOOTHNESS;
  m_rotationSpeed           = CGAL_CLIPPING_PLANE_ROTATION_SPEED;
  m_translationSpeed        = CGAL_CLIPPING_PLANE_TRANSLATION_SPEED;

  m_size = 1.0f;
}

inline 
void Clipping_plane::reset_orientation() 
{
  m_orientation = quatf::Identity(); 
  m_pitch = 0.0f;
  m_targetPitch = 0.0f;
  m_yaw = 0.0f;
  m_targetYaw = 0.0f;

  m_upAxis = vec3f::UnitY();
  m_rightAxis = vec3f::UnitX();
}

inline 
void Clipping_plane::reset_position() 
{
  m_position = vec3f::Zero();
  m_targetPosition = vec3f::Zero();
}

inline 
mat4f Clipping_plane::get_matrix() const
{
  mat4f translation = transform::translation(m_position);
  mat4f rotation = transform::rotation(m_orientation);

  return translation * rotation;
}

inline 
void Clipping_plane::rotation(const float x, const float y) 
{
  if (
    m_constraintAxis == ConstraintAxis::NO_CONSTRAINT ||
    m_constraintAxis == ConstraintAxis::RIGHT_AXIS) 
  {
    m_targetPitch += y;
  }
  
  if (
    m_constraintAxis == ConstraintAxis::NO_CONSTRAINT ||
    m_constraintAxis == ConstraintAxis::UP_AXIS)
  {
    m_targetYaw += x;
  }
}

inline 
void Clipping_plane::translation(const float s)
{
  vec3f forward = (m_orientation * vec3f::UnitZ()).normalized();
  m_targetPosition -= m_size * forward * s * m_translationSpeed; 
}

inline 
void Clipping_plane::translation(const float x, const float y) 
{
  m_targetPosition -= m_size * m_rightAxis * x * m_translationSpeed; 
  m_targetPosition -= m_size * m_upAxis * y * m_translationSpeed; 
}

inline 
void Clipping_plane::translation(const vec3f& direction, const float s) 
{
  m_targetPosition -= m_size * direction.normalized() * s * m_translationSpeed; 
}

inline 
void Clipping_plane::switch_constraint_axis() 
{
  switch(m_constraintAxis)
  {
    case ConstraintAxis::NO_CONSTRAINT: 
      m_constraintAxis = ConstraintAxis::RIGHT_AXIS;
      break;
    case ConstraintAxis::RIGHT_AXIS: 
      m_constraintAxis = ConstraintAxis::UP_AXIS;
      break;
    case ConstraintAxis::UP_AXIS: 
      m_constraintAxis = ConstraintAxis::NO_CONSTRAINT;
      break;
  }
}

#endif // CGAL_CLIPPING_PLANE_H
