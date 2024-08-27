#ifndef CGAL_CAMERA_H
#define CGAL_CAMERA_H

#include <vector>

#include "utils.h"
#include "../bv_settings.h"

class Camera
{
public:
  enum class ConstraintAxis { NO_CONSTRAINT, RIGHT_AXIS, UP_AXIS, FORWARD_AXIS };
  enum class CameraMode { PERSPECTIVE, ORTHOGRAPHIC };
  enum class CameraType { ORBITER, FREE_FLY };

public:
  void update(const float deltaTime);

  void lookat(vec3f const &center, const float size);
  void lookat(vec3f const &pmin, vec3f const &pmax);

  void move(const float z);
  void rotation(const float x, const float y);
  void translation(const float x, const float y);

  void move_up(const float deltaTime);
  void move_down(const float deltaTime);
  void move_right(const float deltaTime);
  void move_left(const float deltaTime);

  mat4f view() const;

  mat4f projection() const;
  mat4f projection(const float width, const float height);

  float znear() const;
  float zfar() const;

  mat4f viewport() const;

  void set_position(const vec3f& position);
  vec3f get_position() const;

  vec3f get_forward() const;
  vec3f get_right() const;
  vec3f get_up() const;

  void reset_all();
  void reset_position();
  void reset_orientation();
  void reset_size();

  void toggle_type(); 

  void switch_constraint_axis();
  void set_constraint_axis(ConstraintAxis axis);

  void increase_fov(float d);
  void disable_smoothness();
  void increase_zoom_smoothness(const float deltaTime);

  void align_to_plane(const vec3f& normal);
  void align_to_nearest_axis();

  std::string get_constraint_axis_str() const;

  bool need_update() const;
  
  inline float get_translation_speed() const { return m_TranslationSpeed; }
  inline float get_rotation_speed() const { return m_RotationSpeed; }

  inline float get_size() const { return m_Size; }
  inline float get_radius() const { return m_Radius; }
  inline vec3f get_center() const { return m_Center; }

  inline float get_fov() const { return m_FOV; }

  inline quatf get_orientation() const { return m_Orientation; }

  inline void set_default_size(float size) { m_DefaultSize = size; }
  inline void set_default_position(const vec3f& position) { m_DefaultPosition = position; set_default_size(position.z()); }
  inline void set_default_orientation(const quatf& orientation) { m_DefaultOrientation = orientation; }

  inline void set_size(float size) { m_Size = size; m_TargetSize = size; }
  inline void set_radius(float radius) { m_Radius = radius; }
  inline void set_center(const vec3f& center) { m_Center = center; }

  inline void set_orthographic() { m_Mode = CameraMode::ORTHOGRAPHIC; }
  inline void set_mode(CameraMode mode) { m_Mode = mode; }
  inline void set_orientation(const quatf& orientation) { m_Orientation = orientation; }
  
  void set_orientation(const vec3f& forward, float upAngle);

  inline void toggle_mode() { m_Mode = (m_Mode == CameraMode::ORTHOGRAPHIC) ? CameraMode::PERSPECTIVE : CameraMode::ORTHOGRAPHIC; }
  inline bool is_orthographic() const { return m_Mode == CameraMode::ORTHOGRAPHIC; }
  inline bool is_orbiter() const { return m_Type == CameraType::ORBITER; }

  inline void increase_rotation_speed(const float deltaTime) { m_RotationSpeed = std::min(m_RotationSpeed+100.f*deltaTime, 360.f); }
  inline void decrease_rotation_speed(const float deltaTime) { m_RotationSpeed = std::max(m_RotationSpeed-100.f*deltaTime, 60.f); }

  inline void increase_translation_speed(const float deltaTime) { m_TranslationSpeed = std::min(m_TranslationSpeed+deltaTime, 20.f); }
  inline void decrease_translation_speed(const float deltaTime) { m_TranslationSpeed = std::max(m_TranslationSpeed-deltaTime, 0.5f); }

  inline void increase_rotation_smoothness(const float deltaTime) { m_RotationSmoothFactor = std::max(m_RotationSmoothFactor-deltaTime, 0.01f); }
  inline void decrease_rotation_smoothness(const float deltaTime) { m_RotationSmoothFactor = std::min(m_RotationSmoothFactor+deltaTime, 1.f); }

  inline void increase_translation_smoothness(const float deltaTime) { m_TranslationSmoothFactor = std::max(m_TranslationSmoothFactor-deltaTime, 0.01f); }
  inline void decrease_translation_smoothness(const float deltaTime) { m_TranslationSmoothFactor = std::min(m_TranslationSmoothFactor+deltaTime, 1.f); }

private: 
  void compute_target_size();

private:
  quatf m_Orientation       { quatf::Identity() }; 
  quatf m_DefaultOrientation  { quatf::Identity() };

  vec3f m_Center { vec3f::Zero() };   

  vec3f m_Position       { vec3f::Zero() }; 
  vec3f m_TargetPosition { vec3f::Zero() };
  vec3f m_DefaultPosition  { vec3f::Zero() }; 

  float m_Size       { CGAL_CAMERA_RADIUS };
  float m_DefaultSize  { CGAL_CAMERA_RADIUS };
  float m_TargetSize { CGAL_CAMERA_RADIUS };
  float m_Radius     { CGAL_CAMERA_RADIUS };

  float m_Width     { 1.0f };  
  float m_Height    { 1.0f }; 
  float m_FOV       { CGAL_CAMERA_FOV };    
  
  float m_RotationSpeed    { CGAL_CAMERA_ROTATION_SPEED }; 
  float m_TranslationSpeed { CGAL_CAMERA_TRANSLATION_SPEED }; 

  CameraType m_Type { CameraType::ORBITER };
  CameraMode m_Mode { CameraMode::PERSPECTIVE };

  ConstraintAxis m_ConstraintAxis { ConstraintAxis::NO_CONSTRAINT };

  float m_Pitch       { 0.0f }; 
  float m_TargetPitch { 0.0f }; 
  float m_Yaw         { 0.0f };
  float m_TargetYaw   { 0.0f }; 
  
  float m_ZoomSmoothFactor        { CGAL_CAMERA_ZOOM_SMOOTHNESS };
  float m_RotationSmoothFactor    { CGAL_CAMERA_ROTATION_SMOOTHNESS };
  float m_TranslationSmoothFactor { CGAL_CAMERA_TRANSLATION_SMOOTHNESS };
};

/********************METHOD IMPLEMENTATIONS********************/

void Camera::update(const float dt) 
{
  if (need_update())
  {
    if (m_Type == CameraType::FREE_FLY) 
    {
      if (m_TargetPitch > 90.f) 
      {
          m_TargetPitch = 90.f;
      }
      else if (m_TargetPitch < -90.f) 
      {
          m_TargetPitch = -90.f;
      }
    }

    float smoothPitch = m_Pitch + m_RotationSmoothFactor * (m_TargetPitch - m_Pitch);   
    float smoothYaw = m_Yaw + m_RotationSmoothFactor * (m_TargetYaw - m_Yaw);   

    float pitchDelta = utils::radians((smoothPitch - m_Pitch) * m_RotationSpeed) * dt;
    float yawDelta = utils::radians((smoothYaw - m_Yaw) * m_RotationSpeed) * dt;

    if (is_orbiter()) 
    {
      if (m_ConstraintAxis == ConstraintAxis::FORWARD_AXIS)
      {
        quatf rollQuaternion(Eigen::AngleAxisf(pitchDelta - yawDelta, -vec3f::UnitZ()));
        m_Orientation = rollQuaternion * m_Orientation;
      }
      else 
      {
        quatf pitchQuaternion(Eigen::AngleAxisf(pitchDelta, vec3f::UnitX()));
        quatf yawQuaternion(Eigen::AngleAxisf(yawDelta, vec3f::UnitY()));
        m_Orientation = pitchQuaternion * yawQuaternion * m_Orientation;
      }
    }
    else // free fly
    {
      quatf pitchQuaternion(Eigen::AngleAxisf(pitchDelta, get_right()));
      quatf yawQuaternion(Eigen::AngleAxisf(yawDelta, vec3f::UnitY())); // lock roll rotation 
      m_Orientation = m_Orientation * pitchQuaternion * yawQuaternion;
    }
    m_Pitch = smoothPitch;
    m_Yaw = smoothYaw;

    m_Position += m_TranslationSmoothFactor * (m_TargetPosition - m_Position);

    m_Size += m_ZoomSmoothFactor * (m_TargetSize - m_Size);
  }
}

void Camera::lookat(vec3f const &center, const float size)
{
  m_Center = center;   // camera focus
  m_Size = size;
  m_DefaultSize = size; // allowing reset
  m_Radius = size;

  compute_target_size();

  m_Orientation = quatf::Identity();
}

void Camera::lookat(vec3f const &pmin, vec3f const &pmax)
{
  lookat(utils::center(pmin, pmax), utils::distance(pmin, pmax));
}

void Camera::move(const float z)
{
  if (is_orbiter()) 
  {
    m_TargetSize -= m_TargetSize * z;
    if (m_TargetSize < 0.001f)
        m_TargetSize = 0.001f;
  }
  else // free fly 
  {
    vec3f forward = get_forward();
    m_TargetPosition += forward * m_Size * z * m_TranslationSpeed;
  }
}

void Camera::rotation(const float x, const float y)
{
  if (
    m_ConstraintAxis == ConstraintAxis::NO_CONSTRAINT ||
    m_ConstraintAxis == ConstraintAxis::RIGHT_AXIS    ||
    m_ConstraintAxis == ConstraintAxis::FORWARD_AXIS) 
  {
    m_TargetPitch += y + y / m_RotationSmoothFactor * .1;
  }
  
  if (
    m_ConstraintAxis == ConstraintAxis::NO_CONSTRAINT ||
    m_ConstraintAxis == ConstraintAxis::UP_AXIS       ||
    m_ConstraintAxis == ConstraintAxis::FORWARD_AXIS)
  {
    m_TargetYaw += x + x / m_RotationSmoothFactor * .1;
  }
}

void Camera::translation(const float x, const float y)
{
  float xspeed = x * m_TranslationSpeed + x / m_TranslationSmoothFactor * .01;
  float yspeed = y * m_TranslationSpeed + y / m_TranslationSmoothFactor * .01;

  if (is_orbiter()) 
  {
    m_TargetPosition.x() += m_Size * xspeed;
    m_TargetPosition.y() += m_Size * yspeed;
  } 
  else // free fly  
  {
    vec3f right = get_right();
    vec3f up = get_up();
    m_TargetPosition -= right * m_Size * xspeed; 
    m_TargetPosition -= up * m_Size * yspeed; 
  }
}

void Camera::move_up(const float dt)
{
  if (is_orbiter()) 
  {
    m_TargetPosition.y() += m_Size * dt * m_TranslationSpeed;
  }
  else // free fly
  {
    vec3f up = get_up();
    m_TargetPosition += up * m_Size * dt * m_TranslationSpeed;
  }
}

void Camera::move_down(const float dt)
{
  move_up(-dt);
}

void Camera::move_right(const float dt)
{
  if (is_orbiter()) 
  {
    m_TargetPosition.x() += m_Size * dt * m_TranslationSpeed;
  }
  else // free fly
  {
    vec3f right = get_right();
    m_TargetPosition += right * m_Size * dt * m_TranslationSpeed;
  }
}

void Camera::move_left(const float dt)
{
  move_right(-dt);
}

mat4f Camera::view() const
{
  mat4f rotation = transform::rotation(m_Orientation);

  if (is_orbiter()) 
  {
    mat4f translation = transform::translation(-m_Position.x(), -m_Position.y(), -m_Size); 
    mat4f model = transform::translation(-m_Center.x(), -m_Center.y(), -m_Center.z()); 
    return translation * rotation * model;
  } 
  else // free fly  
  {
    mat4f translation = transform::translation(-m_Position.x(), -m_Position.y(), -m_Position.z()-m_Size); // translate camera to (0,0,0)
    return rotation * translation;
  }
}

mat4f Camera::projection(const float width, const float height)
{
  m_Width = width;
  m_Height = height;

  return projection();
}

mat4f Camera::projection() const
{ 
  if (m_Mode == CameraMode::ORTHOGRAPHIC) {
    float aspect = m_Width / m_Height;
    float halfWidth = m_Size * aspect * 0.5f;
    float halfHeight = m_Size * 0.5f;
    return utils::ortho(-halfWidth, halfWidth, -halfHeight, halfHeight, znear(), zfar());
  }

  return utils::perspective(utils::radians(m_FOV), m_Width / m_Height, znear(), zfar());
}

float Camera::znear() const
{
  float d;
  if (is_orbiter()) 
  {
    d = utils::distance(m_Center, vec3f(m_Position.x(), m_Position.y(), m_Size));
  }
  else
  {
    d = utils::distance(m_Center, vec3f(m_Position.x(), m_Position.y(), m_Position.z()+m_Size));
  }
  return std::max(0.1f, d - 2 * m_Radius);
}

float Camera::zfar() const
{
  float d;
  if (is_orbiter()) 
  {
    d = utils::distance(m_Center, vec3f(m_Position.x(), m_Position.y(), m_Size));
  }
  else // free fly 
  {
    d = utils::distance(m_Center, vec3f(m_Position.x(), m_Position.y(), m_Position.z()+m_Size));
  }

  return std::max(100.f, d + 2 * m_Radius);
}

mat4f Camera::viewport() const
{
  return transform::viewport(m_Width, m_Height);
}

void Camera::set_position(const vec3f& position) 
{ 
  m_Position = position; 
  m_TargetPosition = position; 
  m_Size = position.z();
  m_TargetSize = position.z();
}

vec3f Camera::get_position() const
{
  if (is_orbiter()) 
  {
    return vec3f(m_Position.x(), m_Position.y(), m_Size);
  }
  else // free fly
  {
    return vec3f(m_Position.x(), m_Position.y(), m_Position.z() + m_Size);
  } 
}

vec3f Camera::get_forward() const
{
  return (m_Orientation.inverse() * -vec3f::UnitZ()).normalized();
}

vec3f Camera::get_right() const
{
  return (m_Orientation.inverse() * vec3f::UnitX()).normalized();
}

vec3f Camera::get_up() const
{
  return (m_Orientation.inverse() * vec3f::UnitY()).normalized();
}

void Camera::reset_all()
{
  reset_position();
  reset_orientation();
  m_FOV                     = CGAL_CAMERA_FOV;
  m_ZoomSmoothFactor        = CGAL_CAMERA_ZOOM_SMOOTHNESS;
  m_RotationSmoothFactor    = CGAL_CAMERA_ROTATION_SMOOTHNESS;
  m_TranslationSmoothFactor = CGAL_CAMERA_TRANSLATION_SMOOTHNESS;
  m_TranslationSpeed        = CGAL_CAMERA_TRANSLATION_SPEED;
  m_RotationSpeed           = CGAL_CAMERA_ROTATION_SPEED;
  reset_size();
}

void Camera::reset_size()
{
  m_Size = m_DefaultSize;
  compute_target_size();
}

void Camera::reset_position()
{
  m_Position = m_DefaultPosition;
  m_TargetPosition = m_DefaultPosition; 
}

void Camera::reset_orientation() 
{ 
  m_Pitch = 0.f; 
  m_Yaw = 0; 
  m_TargetPitch = 0.f; 
  m_TargetYaw = 0; 
  m_Orientation = m_DefaultOrientation; 
}

void Camera::toggle_type() 
{ 
  m_Type = (m_Type == CameraType::ORBITER) ? CameraType::FREE_FLY : CameraType::ORBITER;  
  reset_all(); 
}

void Camera::switch_constraint_axis() 
{
  if (!is_orbiter() && m_ConstraintAxis == ConstraintAxis::UP_AXIS)
  {
    m_ConstraintAxis = ConstraintAxis::NO_CONSTRAINT;
    return;
  }

  switch(m_ConstraintAxis)
  {
    case ConstraintAxis::NO_CONSTRAINT: 
      m_ConstraintAxis = ConstraintAxis::RIGHT_AXIS;
      break;
    case ConstraintAxis::RIGHT_AXIS: 
      m_ConstraintAxis = ConstraintAxis::UP_AXIS;
      break;
    case ConstraintAxis::UP_AXIS: 
      m_ConstraintAxis = ConstraintAxis::FORWARD_AXIS;
      break;
    case ConstraintAxis::FORWARD_AXIS: 
      m_ConstraintAxis = ConstraintAxis::NO_CONSTRAINT;
      break;
  }
}

void Camera::set_constraint_axis(ConstraintAxis axis)
{
  m_ConstraintAxis = axis;
}

void Camera::increase_fov(const float d) 
{ 
  m_FOV += d * 2.f; 
  if (m_FOV > 90.f) m_FOV = 90.;
  if (m_FOV < 45.f) m_FOV = 45.f;

  compute_target_size();
}

void Camera::disable_smoothness()
{
  m_ZoomSmoothFactor = 1.0;
  m_TranslationSmoothFactor = 1.0;
  m_RotationSmoothFactor = 1.0;
}

void Camera::increase_zoom_smoothness(const float s)
{
  m_ZoomSmoothFactor += s * 0.01; 
  if (m_ZoomSmoothFactor > 1.f) m_ZoomSmoothFactor = 1.;
  if (m_ZoomSmoothFactor < 0.01f) m_ZoomSmoothFactor = .01;
}

using NearestAxisResult = std::pair<vec3f, vec3f>;

NearestAxisResult nearest_axis(const vec3f& forward, const vec3f& up) 
{
  std::vector<vec3f> axis = {
    vec3f::UnitX(), -vec3f::UnitX(),
    vec3f::UnitY(), -vec3f::UnitY(),
    vec3f::UnitZ(), -vec3f::UnitZ()
  };

  float maxForwardDot = -1.0f;
  float maxUpDot = -1.0f;
  vec3f nearestForwardAxis = vec3f::Zero();
  vec3f nearestUpAxis = vec3f::Zero();

  for (const auto& a : axis) 
  {
    float dotFw = -forward.dot(a);
    float dotUp = up.dot(a);
    if (dotFw > maxForwardDot) 
    {
      maxForwardDot = dotFw;
      nearestForwardAxis = a; 
    }

    if (dotUp > maxUpDot) 
    {
      maxUpDot = dotUp;
      nearestUpAxis = a; 
    }
  }

  return { nearestForwardAxis, nearestUpAxis };
}

quatf compute_rotation(const vec3f& nearestForwardAxis, const vec3f& nearestUpAxis) 
{
  quatf rotation{};
  if (nearestForwardAxis == vec3f::UnitX()) 
  {
    rotation = quatf(Eigen::AngleAxisf(-M_PI_2, vec3f::UnitY()));
    if (nearestUpAxis == -vec3f::UnitY()) 
    {
      rotation *= quatf(Eigen::AngleAxisf(M_PI, nearestForwardAxis));
    }
    else if (nearestUpAxis == vec3f::UnitZ()) 
    {
      rotation *= quatf(Eigen::AngleAxisf(-M_PI_2, nearestForwardAxis));
    }
    else if (nearestUpAxis == -vec3f::UnitZ()) 
    {
      rotation *= quatf(Eigen::AngleAxisf(M_PI_2, nearestForwardAxis));
    }
  } 
  else if (nearestForwardAxis == -vec3f::UnitX()) 
  {
    rotation = quatf(Eigen::AngleAxisf(M_PI_2, vec3f::UnitY()));
    if (nearestUpAxis == -vec3f::UnitY()) 
    {
      rotation *= quatf(Eigen::AngleAxisf(M_PI, nearestForwardAxis));
    }
    else if (nearestUpAxis == vec3f::UnitZ()) 
    {
      rotation *= quatf(Eigen::AngleAxisf(M_PI_2, nearestForwardAxis));
    }
    else if (nearestUpAxis == -vec3f::UnitZ()) 
    {
      rotation *= quatf(Eigen::AngleAxisf(-M_PI_2, nearestForwardAxis));
    }
  } 
  else if (nearestForwardAxis == vec3f::UnitY()) 
  {
    rotation = quatf(Eigen::AngleAxisf(M_PI_2, vec3f::UnitX()));
    if (nearestUpAxis == vec3f::UnitX()) 
    {
      rotation *= quatf(Eigen::AngleAxisf(M_PI_2, nearestForwardAxis));
    }
    else if (nearestUpAxis == -vec3f::UnitX()) 
    {
      rotation *= quatf(Eigen::AngleAxisf(-M_PI_2, nearestForwardAxis));
    }
    else if (nearestUpAxis == vec3f::UnitZ()) 
    {
      rotation *= quatf(Eigen::AngleAxisf(M_PI, nearestForwardAxis));
    }
  } 
  else if (nearestForwardAxis == -vec3f::UnitY()) 
  {
    rotation = quatf(Eigen::AngleAxisf(-M_PI_2, vec3f::UnitX()));
    if (nearestUpAxis == vec3f::UnitX()) 
    {
      rotation *= quatf(Eigen::AngleAxisf(M_PI_2, nearestForwardAxis));
    }
    else if (nearestUpAxis == -vec3f::UnitX()) 
    {
      rotation *= quatf(Eigen::AngleAxisf(-M_PI_2, nearestForwardAxis));
    }
    else if (nearestUpAxis == -vec3f::UnitZ()) 
    {
      rotation *= quatf(Eigen::AngleAxisf(M_PI, nearestForwardAxis));
    }
  } 
  else if (nearestForwardAxis == vec3f::UnitZ()) 
  {
    rotation = quatf(Eigen::AngleAxisf(0.f, vec3f::UnitY()));
    if (nearestUpAxis == -vec3f::UnitY()) 
    {
      rotation *= quatf(Eigen::AngleAxisf(M_PI, nearestForwardAxis));
    }
    else if (nearestUpAxis == vec3f::UnitX()) 
    {
      rotation *= quatf(Eigen::AngleAxisf(M_PI_2, nearestForwardAxis));
    }
    else if (nearestUpAxis == -vec3f::UnitX()) 
    {
      rotation *= quatf(Eigen::AngleAxisf(-M_PI_2, nearestForwardAxis));
    }
  } 
  else if (nearestForwardAxis == -vec3f::UnitZ()) 
  {
    rotation = quatf(Eigen::AngleAxisf(M_PI, vec3f::UnitY()));
    if (nearestUpAxis == -vec3f::UnitY()) 
    {
      rotation *= quatf(Eigen::AngleAxisf(-M_PI, nearestForwardAxis));
    }
    else if (nearestUpAxis == vec3f::UnitX()) 
    {
      rotation *= quatf(Eigen::AngleAxisf(-M_PI_2, nearestForwardAxis));
    }
    else if (nearestUpAxis == -vec3f::UnitX()) 
    {
      rotation *= quatf(Eigen::AngleAxisf(M_PI_2, nearestForwardAxis));
    }
  }

  return rotation;
}

void Camera::align_to_nearest_axis() 
{
  vec3f forwardDirection = get_forward();
  vec3f upDirection = get_up();

  auto [nearestForwardAxis, nearestUpAxis] = nearest_axis(forwardDirection, upDirection);

  m_Pitch = m_TargetPitch;
  m_Yaw = m_TargetYaw;
  m_Position = m_TargetPosition;

  m_Orientation = compute_rotation(nearestForwardAxis, nearestUpAxis);
}

void Camera::align_to_plane(const vec3f& normal)
{
  vec3f forward = get_forward();
  float dotFF = forward.dot(-normal);   
  float dotFB = forward.dot(normal);  
  
  m_Pitch = m_TargetPitch;
  m_Yaw = m_TargetYaw;
  m_Position = m_TargetPosition;

  if (dotFF > dotFB)
  {
    m_Orientation *= quatf::FromTwoVectors(forward, -normal).inverse();
  }
  else 
  {
    m_Orientation *= quatf::FromTwoVectors(forward,  normal).inverse();
  }
}

void Camera::compute_target_size()
{
  float tanHalfFov = tan(utils::radians(m_FOV) * .5f);
  m_TargetSize = m_Radius / tanHalfFov * .6f;
}

bool Camera::need_update() const
{
  return !utils::equal_float(m_TargetPitch, m_Pitch) 
      || !utils::equal_float(m_TargetYaw, m_Yaw) 
      || !utils::equal_float(m_TargetSize, m_Size) 
      || !utils::equal_vec3f(m_TargetPosition, m_Position)
      ;
}

void Camera::set_orientation(const vec3f& forward, float upAngle) 
{ 
  m_Orientation = quatf::FromTwoVectors(-vec3f::UnitZ(), forward.normalized()).inverse(); 
  m_Orientation *= quatf(Eigen::AngleAxisf(utils::radians(upAngle), get_forward()));
}

std::string Camera::get_constraint_axis_str() const 
{ 
  if (m_ConstraintAxis == ConstraintAxis::UP_AXIS) return "Up";
  if (m_ConstraintAxis == ConstraintAxis::RIGHT_AXIS) return "Right";
  if (m_ConstraintAxis == ConstraintAxis::FORWARD_AXIS) return "Forward";
  return "None"; 
}

#endif // CGAL_CAMERA_H
