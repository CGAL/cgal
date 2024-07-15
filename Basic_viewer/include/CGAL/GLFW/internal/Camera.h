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

  void toggle_type(); 

  void switch_constraint_axis();

  void increase_fov(float d);
  void increase_zoom_smoothness(const float deltaTime);

  void align_to_nearest_axis();

  std::string get_constraint_axis() const;
  
  inline float get_translation_speed() const { return m_translationSpeed; }
  inline float get_rotation_speed() const { return m_rotationSpeed; }

  inline float get_size() const { return m_size; }
  inline float get_radius() const { return m_radius; }
  inline vec3f get_center() const { return m_center; }

  inline quatf get_orientation() const { return m_orientation; }

  inline void set_default_size(float size) { m_defaultSize = size; }
  inline void set_default_position(const vec3f& position) { m_defaultPosition = position; set_default_size(position.z()); }
  inline void set_default_orientation(const quatf& orientation) { m_defaultOrientation = orientation; }

  void set_default_orientation(const vec3f& orientation);

  inline void set_radius(float radius) { m_radius = radius; }
  inline void set_center(const vec3f& center) { m_center = center; }

  inline void set_orthographic() { m_mode = CameraMode::ORTHOGRAPHIC; }
  inline void set_mode(CameraMode mode) { m_mode = mode; }
  inline void set_orientation(const quatf& orientation) { m_orientation = orientation; }
  
  void set_orientation(const vec3f& forward, float upAngle);

  inline void toggle_mode() { m_mode = (m_mode == CameraMode::ORTHOGRAPHIC) ? CameraMode::PERSPECTIVE : CameraMode::ORTHOGRAPHIC; }
  inline bool is_orthographic() const { return m_mode == CameraMode::ORTHOGRAPHIC; }
  inline bool is_orbiter() const { return m_type == CameraType::ORBITER; }

  inline void increase_rotation_speed(const float deltaTime) { m_rotationSpeed = std::min(m_rotationSpeed+100.f*deltaTime, 360.f); }
  inline void decrease_rotation_speed(const float deltaTime) { m_rotationSpeed = std::max(m_rotationSpeed-100.f*deltaTime, 60.f); }

  inline void increase_translation_speed(const float deltaTime) { m_translationSpeed = std::min(m_translationSpeed+deltaTime, 20.f); }
  inline void decrease_translation_speed(const float deltaTime) { m_translationSpeed = std::max(m_translationSpeed-deltaTime, 0.5f); }

  inline void increase_rotation_smoothness(const float deltaTime) { m_rotationSmoothFactor = std::max(m_rotationSmoothFactor-deltaTime, 0.01f); }
  inline void decrease_rotation_smoothness(const float deltaTime) { m_rotationSmoothFactor = std::min(m_rotationSmoothFactor+deltaTime, 1.f); }

  inline void increase_translation_smoothness(const float deltaTime) { m_translationSmoothFactor = std::max(m_translationSmoothFactor-deltaTime, 0.01f); }
  inline void decrease_translation_smoothness(const float deltaTime) { m_translationSmoothFactor = std::min(m_translationSmoothFactor+deltaTime, 1.f); }

private: 
  void compute_target_size();

private:
  quatf m_orientation       { quatf::Identity() }; 
  quatf m_defaultOrientation  { quatf::Identity() };

  vec3f m_center { vec3f::Zero() };   

  vec3f m_position       { vec3f::Zero() }; 
  vec3f m_targetPosition { vec3f::Zero() };
  vec3f m_defaultPosition  { vec3f::Zero() }; 

  float m_size       { CGAL_CAMERA_RADIUS };
  float m_defaultSize  { CGAL_CAMERA_RADIUS };
  float m_targetSize { CGAL_CAMERA_RADIUS };
  float m_radius     { CGAL_CAMERA_RADIUS };

  float m_width     { 1.0f };  
  float m_height    { 1.0f }; 
  float m_fov       { CGAL_CAMERA_FOV };    
  
  float m_rotationSpeed    { CGAL_CAMERA_ROTATION_SPEED }; 
  float m_translationSpeed { CGAL_CAMERA_TRANSLATION_SPEED }; 

  CameraType m_type { CameraType::ORBITER };
  CameraMode m_mode { CameraMode::PERSPECTIVE };

  ConstraintAxis m_constraintAxis { ConstraintAxis::NO_CONSTRAINT };

  float m_pitch       { 0.0f }; 
  float m_targetPitch { 0.0f }; 
  float m_yaw         { 0.0f };
  float m_targetYaw   { 0.0f }; 
  
  float m_zoomSmoothFactor        { CGAL_CAMERA_ZOOM_SMOOTHNESS };
  float m_rotationSmoothFactor    { CGAL_CAMERA_ROTATION_SMOOTHNESS };
  float m_translationSmoothFactor { CGAL_CAMERA_TRANSLATION_SMOOTHNESS };
};

/********************METHOD DEFINITIONS********************/

inline 
void Camera::update(const float dt) 
{
  if (m_type == CameraType::FREE_FLY) {
    if (m_targetPitch > 90.f) 
    {
        m_targetPitch = 90.f;
    }
    else if (m_targetPitch < -90.f) 
    {
        m_targetPitch = -90.f;
    }
  }

  float smoothPitch = m_pitch + m_rotationSmoothFactor * (m_targetPitch - m_pitch);   
  float smoothYaw = m_yaw + m_rotationSmoothFactor * (m_targetYaw - m_yaw);   

  float pitchDelta = utils::radians((smoothPitch - m_pitch) * m_rotationSpeed) * dt;
  float yawDelta = utils::radians((smoothYaw - m_yaw) * m_rotationSpeed) * dt;

  if (is_orbiter()) 
  {
    if (m_constraintAxis == ConstraintAxis::FORWARD_AXIS)
    {
      quatf rollQuaternion(Eigen::AngleAxisf(pitchDelta - yawDelta, -vec3f::UnitZ()));
      m_orientation = rollQuaternion * m_orientation;
    }
    else 
    {
      quatf pitchQuaternion(Eigen::AngleAxisf(pitchDelta, vec3f::UnitX()));
      quatf yawQuaternion(Eigen::AngleAxisf(yawDelta, vec3f::UnitY()));
      m_orientation = pitchQuaternion * yawQuaternion * m_orientation;
    }

  }
  else // free fly
  {
    quatf pitchQuaternion(Eigen::AngleAxisf(pitchDelta, get_right()));
    quatf yawQuaternion(Eigen::AngleAxisf(yawDelta, get_up()));
    m_orientation = m_orientation * pitchQuaternion * yawQuaternion;
  }
  m_pitch = smoothPitch;
  m_yaw = smoothYaw;

  m_position += m_translationSmoothFactor * (m_targetPosition - m_position);

  m_size += m_zoomSmoothFactor * (m_targetSize - m_size);
}

inline 
void Camera::lookat(vec3f const &center, const float size)
{
  m_center = center;   // camera focus
  m_size = size;
  m_defaultSize = size; // allowing reset
  m_radius = size;

  compute_target_size();

  m_orientation = quatf::Identity();
}

inline 
void Camera::lookat(vec3f const &pmin, vec3f const &pmax)
{
  lookat(utils::center(pmin, pmax), utils::distance(pmin, pmax));
}

inline 
void Camera::move(const float z)
{
  if (is_orbiter()) 
  {
    m_targetSize -= m_targetSize * z;
    if (m_targetSize < 0.001f)
        m_targetSize = 0.001f;
  }
  else // free fly 
  {
    vec3f forward = get_forward();
    m_targetPosition += forward * m_size * z * m_translationSpeed;
  }
}

inline 
void Camera::rotation(const float x, const float y)
{
  if (
    m_constraintAxis == ConstraintAxis::NO_CONSTRAINT ||
    m_constraintAxis == ConstraintAxis::RIGHT_AXIS    ||
    m_constraintAxis == ConstraintAxis::FORWARD_AXIS) 
  {
    m_targetPitch += y + y / m_rotationSmoothFactor * .1;
  }
  
  if (
    m_constraintAxis == ConstraintAxis::NO_CONSTRAINT ||
    m_constraintAxis == ConstraintAxis::UP_AXIS       ||
    m_constraintAxis == ConstraintAxis::FORWARD_AXIS)
  {
    m_targetYaw += x + x / m_rotationSmoothFactor * .1;
  }
}

inline 
void Camera::translation(const float x, const float y)
{
  float xspeed = x * m_translationSpeed + x / m_translationSmoothFactor * .01;
  float yspeed = y * m_translationSpeed + y / m_translationSmoothFactor * .01;

  if (is_orbiter()) 
  {
    m_targetPosition.x() += m_size * xspeed;
    m_targetPosition.y() += m_size * yspeed;
  } 
  else // free fly  
  {
    vec3f right = get_right();
    vec3f up = get_up();
    m_targetPosition -= right * m_size * xspeed; 
    m_targetPosition -= up * m_size * yspeed; 
  }
}

inline 
void Camera::move_up(const float dt)
{
  if (is_orbiter()) 
  {
    m_targetPosition.y() += m_size * dt * m_translationSpeed;
  }
  else // free fly
  {
    vec3f up = get_up();
    m_targetPosition += up * m_size * dt * m_translationSpeed;
  }
}

inline 
void Camera::move_down(const float dt)
{
  move_up(-dt);
}

inline 
void Camera::move_right(const float dt)
{
  if (is_orbiter()) 
  {
    m_targetPosition.x() += m_size * dt * m_translationSpeed;
  }
  else // free fly
  {
    vec3f right = get_right();
    m_targetPosition += right * m_size * dt * m_translationSpeed;
  }
}

inline 
void Camera::move_left(const float dt)
{
  move_right(-dt);
}

inline 
mat4f Camera::view() const
{
  mat4f rotation = transform::rotation(m_orientation);

  if (is_orbiter()) 
  {
    mat4f translation = transform::translation(-m_position.x(), -m_position.y(), -m_size); 
    mat4f model = transform::translation(-m_center.x(), -m_center.y(), -m_center.z()); 
    return translation * rotation * model;
  } 
  else // free fly  
  {
    mat4f translation = transform::translation(-m_position.x(), -m_position.y(), -m_position.z()-m_size); // translate camera to (0,0,0)
    return rotation * translation;
  }
}

inline 
mat4f Camera::projection(const float width, const float height)
{
  m_width = width;
  m_height = height;

  return projection();
}

inline 
mat4f Camera::projection() const
{ 
  if (m_mode == CameraMode::ORTHOGRAPHIC) {
    float aspect = m_width / m_height;
    float halfWidth = m_size * aspect * 0.5f;
    float halfHeight = m_size * 0.5f;
    return utils::ortho(-halfWidth, halfWidth, -halfHeight, halfHeight, znear(), zfar());
  }

  return utils::perspective(utils::radians(m_fov), m_width / m_height, znear(), zfar());
}

inline 
float Camera::znear() const
{
  float d;
  if (is_orbiter()) 
  {
    d = utils::distance(m_center, vec3f(m_position.x(), m_position.y(), m_size));
  }
  else
  {
    d = utils::distance(m_center, vec3f(m_position.x(), m_position.y(), m_position.z()+m_size));
  }
  return std::max(0.1f, d - 2 * m_radius);
}

inline 
float Camera::zfar() const
{
  float d;
  if (is_orbiter()) 
  {
    d = utils::distance(m_center, vec3f(m_position.x(), m_position.y(), m_size));
  }
  else // free fly 
  {
    d = utils::distance(m_center, vec3f(m_position.x(), m_position.y(), m_position.z()+m_size));
  }

  return std::max(100.f, d + 2 * m_radius);
}

inline 
mat4f Camera::viewport() const
{
  return transform::viewport(m_width, m_height);
}

inline 
void Camera::set_position(const vec3f& position) 
{ 
  m_position = position; 
  m_targetPosition = position; 
  m_size = position.z();
  m_targetSize = position.z();
}

inline 
vec3f Camera::get_position() const
{
  if (is_orbiter()) 
  {
    return vec3f(m_position.x(), m_position.y(), m_size);
  }
  else // free fly
  {
    return vec3f(m_position.x(), m_position.y(), m_position.z() + m_size);
  } 
}

inline 
vec3f Camera::get_forward() const
{
  vec3f forward;
  forward = m_orientation.inverse() * -vec3f::UnitZ();

  return forward.normalized();
}

inline 
vec3f Camera::get_right() const
{
  vec3f right;
  right = m_orientation.inverse() * vec3f::UnitX();

  return right.normalized();
}

inline 
vec3f Camera::get_up() const
{
  vec3f up;
  up = m_orientation.inverse() * vec3f::UnitY();

  return up.normalized();
}

inline 
void Camera::reset_all()
{
  reset_position();
  reset_orientation();
  m_size = m_defaultSize;
  m_fov                     = CGAL_CAMERA_FOV;
  m_zoomSmoothFactor        = CGAL_CAMERA_ZOOM_SMOOTHNESS;
  m_rotationSmoothFactor    = CGAL_CAMERA_ROTATION_SMOOTHNESS;
  m_translationSmoothFactor = CGAL_CAMERA_TRANSLATION_SMOOTHNESS;
  m_translationSpeed        = CGAL_CAMERA_TRANSLATION_SPEED;
  m_rotationSpeed           = CGAL_CAMERA_ROTATION_SPEED;

  compute_target_size();
}

inline 
void Camera::reset_position()
{
  m_position = m_defaultPosition;
  m_targetPosition = m_defaultPosition; 
}

inline 
void Camera::reset_orientation() 
{ 
  m_pitch = 0.f; 
  m_yaw = 0; 
  m_targetPitch = 0.f; 
  m_targetYaw = 0; 
  m_orientation = m_defaultOrientation; 
}

inline 
void Camera::toggle_type() 
{ 
  m_type = (m_type == CameraType::ORBITER) ? CameraType::FREE_FLY : CameraType::ORBITER;  
  reset_all(); 
}

inline 
void Camera::switch_constraint_axis() 
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
      m_constraintAxis = ConstraintAxis::FORWARD_AXIS;
      break;
    case ConstraintAxis::FORWARD_AXIS: 
      m_constraintAxis = ConstraintAxis::NO_CONSTRAINT;
      break;
  }
}

inline 
void Camera::increase_fov(const float d) 
{ 
  m_fov += d * 2.f; 
  if (m_fov > 90.f) m_fov = 90.;
  if (m_fov < 45.f) m_fov = 45.f;

  compute_target_size();
}

inline 
void Camera::increase_zoom_smoothness(const float s)
{
  m_zoomSmoothFactor += s * 0.01; 
  if (m_zoomSmoothFactor > 1.f) m_zoomSmoothFactor = 1.;
  if (m_zoomSmoothFactor < 0.01f) m_zoomSmoothFactor = .01;
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

inline 
void Camera::align_to_nearest_axis() 
{
  vec3f forwardDirection = get_forward();
  vec3f upDirection = get_up();

  auto [nearestForwardAxis, nearestUpAxis] = nearest_axis(forwardDirection, upDirection);

  m_pitch = m_targetPitch;
  m_yaw = m_targetYaw;
  m_position = m_targetPosition;

  m_orientation = compute_rotation(nearestForwardAxis, nearestUpAxis);
}

inline
void Camera::compute_target_size()
{
  float tanHalfFov = tan(utils::radians(m_fov) * .5f);
  m_targetSize = m_radius / tanHalfFov * .6f;
}

inline 
void Camera::set_orientation(const vec3f& forward, float upAngle) 
{ 
  m_orientation = quatf::FromTwoVectors(-vec3f::UnitZ(), forward.normalized()).inverse(); 
  m_orientation *= quatf(Eigen::AngleAxisf(utils::radians(upAngle), get_forward()));
}

inline 
std::string Camera::get_constraint_axis() const 
{ 
  if (m_constraintAxis == ConstraintAxis::UP_AXIS) return "Up";
  if (m_constraintAxis == ConstraintAxis::RIGHT_AXIS) return "Right";
  if (m_constraintAxis == ConstraintAxis::FORWARD_AXIS) return "Forward";
  return "None"; 
}

#endif // CGAL_CAMERA_H
