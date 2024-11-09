#ifndef CGAL_GLFW_INTERNAL_CLIPPING_PLANE_H
#define CGAL_GLFW_INTERNAL_CLIPPING_PLANE_H

#include "utils.h"
#include "Line_renderer.h"
#include "../bv_settings.h"

class Clipping_plane : public Line_renderer
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

  vec3f get_normal() const;

  bool need_update() const;

  inline float get_transparency() const { return m_Transparency; }
  inline float get_rotation_speed() const { return m_RotationSpeed; }
  inline float get_translation_speed() const { return m_TranslationSpeed; }
  inline std::string get_constraint_axis_str() const { return m_ConstraintAxis == ConstraintAxis::NO_CONSTRAINT ? "None" : (m_ConstraintAxis == ConstraintAxis::RIGHT_AXIS ? "Right" : "Up"); }

  void set_orientation(const vec3f& normal);

  inline void set_size(const float size) { m_Size = size; }
  inline void set_up_axis(const vec3f& upAxis) { m_UpAxis = upAxis; }
  inline void set_right_axis(const vec3f& rightAxis) { m_RightAxis = rightAxis; }

  inline void increase_rotation_speed(const float deltaTime) { m_RotationSpeed = std::min(m_RotationSpeed+100.f*deltaTime, 360.f); }
  inline void decrease_rotation_speed(const float deltaTime) { m_RotationSpeed = std::max(m_RotationSpeed-100.f*deltaTime, 60.f); }

  inline void increase_translation_speed(const float deltaTime) { m_TranslationSpeed = std::min(m_TranslationSpeed+deltaTime, 20.f); }
  inline void decrease_translation_speed(const float deltaTime) { m_TranslationSpeed = std::max(m_TranslationSpeed-deltaTime, 0.5f); }

  inline void increase_rotation_smoothness(const float deltaTime) { m_RotationSmoothFactor = std::max(m_RotationSmoothFactor-deltaTime, 0.01f); }
  inline void decrease_rotation_smoothness(const float deltaTime) { m_RotationSmoothFactor = std::min(m_RotationSmoothFactor+deltaTime, 1.f); }

  inline void increase_translation_smoothness(const float deltaTime) { m_TranslationSmoothFactor = std::max(m_TranslationSmoothFactor-deltaTime, 0.01f); }
  inline void decrease_translation_smoothness(const float deltaTime) { m_TranslationSmoothFactor = std::min(m_TranslationSmoothFactor+deltaTime, 1.f); }

  void align_to_direction(const vec3f& direction);

private:
  quatf m_Orientation { quatf::Identity() };

  vec3f m_Position       { vec3f::Zero() };
  vec3f m_TargetPosition { vec3f::Zero() };

  vec3f m_UpAxis    { vec3f::UnitY() };
  vec3f m_RightAxis { vec3f::UnitX() };

  float m_Pitch       { 0.0f };
  float m_TargetPitch { 0.0f };
  float m_Yaw         { 0.0f };
  float m_TargetYaw   { 0.0f };

  float m_RotationSpeed    { CGAL_CLIPPING_PLANE_ROTATION_SPEED };
  float m_TranslationSpeed { CGAL_CLIPPING_PLANE_TRANSLATION_SPEED };

  float m_RotationSmoothFactor    { CGAL_CLIPPING_PLANE_ROTATION_SMOOTHNESS };
  float m_TranslationSmoothFactor { CGAL_CLIPPING_PLANE_TRANSLATION_SMOOTHNESS };

  float m_Size { 1.0f }; 

  float m_Transparency { CGAL_CLIPPING_PLANE_RENDERING_TRANSPARENCY }; // to what extent the transparent part should be rendered;

  ConstraintAxis m_ConstraintAxis { ConstraintAxis::NO_CONSTRAINT };
};

/********************METHOD IMPLEMENTATIONS********************/
 
void Clipping_plane::update(const float deltaTime)
{
  if (need_update())
  {
    float smoothPitch = m_Pitch + m_RotationSmoothFactor * (m_TargetPitch - m_Pitch);   
    float smoothYaw = m_Yaw + m_RotationSmoothFactor * (m_TargetYaw - m_Yaw);   

    float pitchDelta = utils::radians((smoothPitch - m_Pitch) * m_RotationSpeed) * deltaTime;
    float yawDelta = utils::radians((smoothYaw - m_Yaw) * m_RotationSpeed) * deltaTime;

    quatf pitchQuaternion(Eigen::AngleAxisf(pitchDelta, m_RightAxis));
    quatf yawQuaternion(Eigen::AngleAxisf(yawDelta, m_UpAxis));
    m_Orientation = pitchQuaternion * yawQuaternion * m_Orientation;

    m_Pitch = smoothPitch;
    m_Yaw = smoothYaw;

    m_Position += m_TranslationSmoothFactor * (m_TargetPosition - m_Position);
  }
}
 
void Clipping_plane::reset_all()
{
  reset_position();
  reset_orientation();
  
  m_RotationSmoothFactor    = CGAL_CLIPPING_PLANE_ROTATION_SMOOTHNESS;
  m_TranslationSmoothFactor = CGAL_CLIPPING_PLANE_TRANSLATION_SMOOTHNESS;
  m_RotationSpeed           = CGAL_CLIPPING_PLANE_ROTATION_SPEED;
  m_TranslationSpeed        = CGAL_CLIPPING_PLANE_TRANSLATION_SPEED;

  m_Size = 1.0f;
}
 
void Clipping_plane::reset_orientation() 
{
  m_Orientation = quatf::Identity(); 
  m_Pitch = 0.0f;
  m_TargetPitch = 0.0f;
  m_Yaw = 0.0f;
  m_TargetYaw = 0.0f;

  m_UpAxis = vec3f::UnitY();
  m_RightAxis = vec3f::UnitX();
}
 
void Clipping_plane::reset_position() 
{
  m_Position = vec3f::Zero();
  m_TargetPosition = vec3f::Zero();
}
 vec3f Clipping_plane::get_normal() const
{
  return (m_Orientation * vec3f::UnitZ()).normalized();
}
 
mat4f Clipping_plane::get_matrix() const
{
  mat4f translation = transform::translation(m_Position);
  mat4f rotation = transform::rotation(m_Orientation);

  return translation * rotation;
}
 
void Clipping_plane::rotation(const float x, const float y) 
{
  if (
    m_ConstraintAxis == ConstraintAxis::NO_CONSTRAINT ||
    m_ConstraintAxis == ConstraintAxis::RIGHT_AXIS) 
  {
    m_TargetPitch += y;
  }
  
  if (
    m_ConstraintAxis == ConstraintAxis::NO_CONSTRAINT ||
    m_ConstraintAxis == ConstraintAxis::UP_AXIS)
  {
    m_TargetYaw += x;
  }
}
 
void Clipping_plane::translation(const float s)
{
  vec3f normal = get_normal();
  m_TargetPosition -= m_Size * normal * s * m_TranslationSpeed; 
}
 
void Clipping_plane::translation(const float x, const float y) 
{
  m_TargetPosition -= m_Size * m_RightAxis * x * m_TranslationSpeed; 
  m_TargetPosition -= m_Size * m_UpAxis * y * m_TranslationSpeed; 
}
 
void Clipping_plane::translation(const vec3f& direction, const float s) 
{
  m_TargetPosition -= m_Size * direction.normalized() * s * m_TranslationSpeed; 
}
 
void Clipping_plane::switch_constraint_axis() 
{
  switch(m_ConstraintAxis)
  {
    case ConstraintAxis::NO_CONSTRAINT: 
      m_ConstraintAxis = ConstraintAxis::RIGHT_AXIS;
      break;
    case ConstraintAxis::RIGHT_AXIS: 
      m_ConstraintAxis = ConstraintAxis::UP_AXIS;
      break;
    case ConstraintAxis::UP_AXIS: 
      m_ConstraintAxis = ConstraintAxis::NO_CONSTRAINT;
      break;
  }
}
 
bool Clipping_plane::need_update() const
{
  return !utils::equal_float(m_TargetPitch, m_Pitch) 
      || !utils::equal_float(m_TargetYaw, m_Yaw) 
      || !utils::equal_vec3f(m_TargetPosition, m_Position)
      ; 
}
 
void Clipping_plane::set_orientation(const vec3f& normal)
{
  m_Orientation = quatf::FromTwoVectors(vec3f::UnitZ(), normal);
}

void Clipping_plane::align_to_direction(const vec3f& direction)
{
  vec3f normal = get_normal();
  float dotFF = normal.dot(-direction);   
  float dotFB = normal.dot(direction);  
  
  m_Pitch = m_TargetPitch;
  m_Yaw = m_TargetYaw;
  m_Position = m_TargetPosition;

  if (dotFF > dotFB)
  {
    m_Orientation = quatf::FromTwoVectors(normal, -direction) * m_Orientation;
  }
  else 
  {
    m_Orientation = quatf::FromTwoVectors(normal,  direction) * m_Orientation;
  }
}

#endif // CGAL_GLFW_INTERNAL_CLIPPING_PLANE_H
