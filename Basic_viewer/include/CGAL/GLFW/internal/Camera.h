#ifndef CGAL_CAMERA_H
#define CGAL_CAMERA_H

#include "utils.h"

class Camera
{
public:
  enum class Mode
  {
    PERSPECTIVE,
    ORTHOGRAPHIC,
  };

  enum class Type
  {
    ORBITER,
    FREE_FLY
  };

public:
  Camera();

  void lookat(vec3f const &center, const float size);
  void lookat(vec3f const &pmin, vec3f const &pmax);

  void update(const float dt);

  void translation(const float x, const float y);
  void rotation(const float x, const float y);
  void move(const float z);

  void move_up(const float dt);
  void move_down(const float dt);
  void move_right(const float dt);
  void move_left(const float dt);

  mat4f view() const;

  mat4f projection() const;
  mat4f projection(const float width, const float height);

  float znear() const;
  float zfar() const;

  mat4f viewport() const;

  mat4f image_toWorld() const;
  void frame(const float z, vec3f &dO, vec3f &dx, vec3f &dy) const;

  vec3f position() const;
  vec3f forward_direction() const;
  vec3f right_direction() const;
  vec3f up_direction() const;

  void set_fov(float d);

  inline float get_size() const { return m_size; }

  void reset_position();
  void reset_rotation();
  void reset_all();

  void align_to_nearest_axis();

  inline float get_radius() const { return m_radius; }
  inline vec3f get_center() const { return m_center; }

  inline quatf get_orientation() const { return m_orientation; }

  void toggle_fly(); 
  void set_zoom_smoothness(const float s);

  inline void mode(Mode mode) { m_mode = mode; }
  inline void toggle_mode() { m_mode = (m_mode == Mode::ORTHOGRAPHIC) ? Mode::PERSPECTIVE : Mode::ORTHOGRAPHIC; }
  inline bool is_orthographic() const { return m_mode == Mode::ORTHOGRAPHIC; }
  inline bool is_orbiter() const { return m_type == Type::ORBITER; }

  inline void inc_rspeed() { m_rotationSpeed = std::min(m_rotationSpeed+10.f, 500.f); }
  inline void dec_rspeed() { m_rotationSpeed = std::max(m_rotationSpeed-10.f, 50.f); }

  inline void inc_tspeed() { m_translationSpeed = std::min(m_translationSpeed+.1f, 20.f); }
  inline void dec_tspeed() { m_translationSpeed = std::max(m_translationSpeed-.1f, 0.5f); }

  inline void inc_rotation_smoothness() { m_rotationSmoothFactor = std::max(m_rotationSmoothFactor-.01f, 0.01f); }
  inline void dec_rotation_smoothness() { m_rotationSmoothFactor = std::min(m_rotationSmoothFactor+.01f, 1.f); }

  inline void inc_translation_smoothness() { m_translationSmoothFactor = std::max(m_translationSmoothFactor-.01f, 0.01f); }
  inline void dec_translation_smoothness() { m_translationSmoothFactor = std::min(m_translationSmoothFactor+.01f, 1.f); }

private:
  vec3f m_center;   // centre de l'objet observé
  vec3f m_position; // position de la caméra (translations appliquées)

  quatf m_orientation; // orientation de la caméra (rotations appliquées)

  float m_size;
  float m_constSize;
  float m_radius;

  float m_width;  // largeur de l'image
  float m_height; // hauteur de l'image
  float m_fov;    // FOV pour la projection
  
  float m_translationSpeed; 
  float m_rotationSpeed; 

  Type m_type;
  Mode m_mode;

  /* SMOOTHNESS PARAMS */

  float m_targetSize;

  float m_pitch; 
  float m_yaw; 
  
  float m_targetPitch; 
  float m_targetYaw; 

  float m_rotationSmoothFactor;
  float m_translationSmoothFactor;
  float m_zoomSmoothFactor;

  vec3f m_targetPosition;
};

/********************DEFINITIONS********************/

Camera::Camera() : 
  m_type(Type::ORBITER),
  m_mode(Mode::PERSPECTIVE), 
  m_center(), 
  m_position(vec3f::Zero()), 
  m_targetPosition(vec3f::Zero()), 
  m_orientation(quatf::Identity()), 
  m_translationSpeed(1.f), m_rotationSpeed(100.f), 
  m_constSize(5.f), 
  m_radius(5.f),
  m_width(1.f), m_height(1.f),
  m_pitch(0.f),
  m_targetPitch(0.f),
  m_yaw(0.f),
  m_targetYaw(0.f),
  m_size(5.f),
  m_targetSize(5.f),
  m_fov(45.f),
  m_zoomSmoothFactor(0.1f),
  m_rotationSmoothFactor(0.3f),
  m_translationSmoothFactor(0.3f)
{

}

inline 
void Camera::lookat(vec3f const &center, const float size)
{
  m_center = center;   // camera focus
  m_size = size;
  m_targetSize = size;
  m_constSize = size; // allowing reset
  m_radius = size;

  m_orientation = quatf::Identity();
}

inline 
void Camera::lookat(vec3f const &pmin, vec3f const &pmax)
{
  lookat(center(pmin, pmax), distance(pmin, pmax));
}

inline 
void Camera::toggle_fly() 
{ 
  m_type = (m_type == Type::ORBITER) ? Type::FREE_FLY : Type::ORBITER;  
  reset_all(); 
}

inline 
void Camera::update(const float dt) 
{
  if (m_type == Type::FREE_FLY) {
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

  float pitchDelta = radians((smoothPitch - m_pitch) * m_rotationSpeed) * dt;
  float yawDelta = radians((smoothYaw - m_yaw) * m_rotationSpeed) * dt;

  if (is_orbiter()) 
  {
    quatf pitchQuaternion(Eigen::AngleAxisf(pitchDelta, vec3f::UnitX()));
    quatf yawQuaternion(Eigen::AngleAxisf(yawDelta, vec3f::UnitY()));
    m_orientation = pitchQuaternion * yawQuaternion * m_orientation;
  }
  else // free fly
  {
    quatf pitchQuaternion(Eigen::AngleAxisf(pitchDelta, right_direction()));
    quatf yawQuaternion(Eigen::AngleAxisf(yawDelta, vec3f::UnitY()));
    m_orientation = m_orientation * pitchQuaternion * yawQuaternion;
  }
  m_pitch = smoothPitch;
  m_yaw = smoothYaw;

  m_position += m_translationSmoothFactor * (m_targetPosition - m_position);

  m_size += m_zoomSmoothFactor * (m_targetSize - m_size);
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
    vec3f right = right_direction();
    vec3f up = up_direction();
    m_targetPosition -= right * m_size * xspeed; 
    m_targetPosition -= up * m_size * yspeed; 
  }
}

inline 
void Camera::rotation(const float x, const float y)
{
  m_targetPitch += y + y / m_rotationSmoothFactor * .1;
  m_targetYaw += x + x / m_rotationSmoothFactor * .1;
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
    vec3f forward = forward_direction();
    m_targetPosition -= forward * m_size * z * m_translationSpeed;
  }
}

inline 
void Camera::move_left(const float dt)
{
  move_right(-dt);
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
    vec3f right = right_direction();
    m_targetPosition += right * m_size * dt * m_translationSpeed;
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
    vec3f up = up_direction();
    m_targetPosition += up * m_size * dt * m_translationSpeed;
  }
}

inline 
void Camera::move_down(const float dt)
{
  move_up(-dt);
}

inline 
mat4f Camera::view() const
{
  mat3f rotation3x3 = m_orientation.toRotationMatrix();
  mat4f rotation = mat4f::Identity();
  rotation.block<3, 3>(0, 0) = rotation3x3;

  if (is_orbiter()) 
  {
    mat4f camera = transform::translation(-m_position.x(), -m_position.y(), -m_size); 
    mat4f model = transform::translation(-m_center.x(), -m_center.y(), -m_center.z()); 
    return  camera * rotation * model;
  } 
  else // free fly  
  {
    mat4f translation = transform::translation(-m_position.x(), -m_position.y(), -m_position.z()-m_size); // translate camera to (0,0,0)
    return rotation * translation;
  }
}

inline 
void Camera::set_fov(const float d) 
{ 
  m_fov += d * 2.f; 
  if (m_fov > 120.f) m_fov = 120.;
  if (m_fov < 15.f) m_fov = 15.f;

  float aspectRatio = m_width / m_height;
  float tanHalfFov = tan(radians(m_fov) * .5f);
  m_targetSize = m_radius / tanHalfFov * .6f;
}

inline 
void Camera::set_zoom_smoothness(const float s)
{
  m_zoomSmoothFactor += s * 0.01; 
  if (m_zoomSmoothFactor > 1.f) m_zoomSmoothFactor = 1.;
  if (m_zoomSmoothFactor < 0.01f) m_zoomSmoothFactor = .01;
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
  if (m_mode == Mode::ORTHOGRAPHIC) {
    float aspect = m_width / m_height;
    float halfWidth = m_size * aspect * 0.5f;
    float halfHeight = m_size * 0.5f;
    return ortho(-halfWidth, halfWidth, -halfHeight, halfHeight, znear(), zfar());
  }

  return PERSPECTIVE(radians(m_fov), m_width / m_height, znear(), zfar());
}

inline 
float Camera::znear() const
{
  float d;
  if (is_orbiter()) 
  {
    d = distance(m_center, vec3f(m_position.x(), m_position.y(), m_size));
  }
  else
  {
    d = distance(m_center, vec3f(m_position.x(), m_position.y(), m_position.z()+m_size));
  }
  return std::max(0.1f, d - 2 * m_radius);
}

inline 
float Camera::zfar() const
{
  float d;
  if (is_orbiter()) 
  {
    d = distance(m_center, vec3f(m_position.x(), m_position.y(), m_size));
  }
  else
  {
    d = distance(m_center, vec3f(m_position.x(), m_position.y(), m_position.z()+m_size));
  }

  return std::max(100.f, d + 2 * m_radius);
}

inline 
mat4f Camera::viewport() const
{
  return transform::viewport(m_width, m_height);
}

inline 
vec3f Camera::position() const
{
  mat4f v = view();
  mat4f vi = v.inverse();

  vec3f o(0.f, 0.f, 0.f);
  return mult_vec_mat(o, vi);
}

inline 
vec3f Camera::forward_direction() const
{
  vec3f forward;
  if (is_orbiter()) 
  {
    forward.x() = m_center.x() - position().x();
    forward.y() = m_center.y() - position().y();
    forward.z() = m_center.z() - position().z();
  }
  else 
  {
    forward = m_orientation.inverse() * vec3f::UnitZ();
  }

  return forward.normalized();
}

inline 
vec3f Camera::right_direction() const
{
  vec3f right;
  right = m_orientation.inverse() * vec3f::UnitX();

  return right.normalized();
}

inline 
vec3f Camera::up_direction() const
{
  vec3f up;
  up = m_orientation.inverse() * vec3f::UnitY();

  return up.normalized();
}

inline 
void Camera::reset_rotation() 
{ 
  m_pitch = 0.f; 
  m_yaw = 0; 
  m_targetPitch = 0.f; 
  m_targetYaw = 0; 
  m_orientation = quatf::Identity(); 
}


inline void Camera::reset_position()
{
  m_position = vec3f::Zero();
  m_targetPosition = vec3f::Zero(); 
}


inline 
void Camera::reset_all()
{
  reset_position();
  reset_rotation();
  m_size = m_constSize;
  m_targetSize = m_constSize;
  m_fov = 45.f;
  m_zoomSmoothFactor = .1f;
  m_rotationSmoothFactor = .1f;
  m_translationSmoothFactor = .1f;
}

inline 
mat4f Camera::image_toWorld() const 
{
  mat4f v = view();
  mat4f p = projection();
  mat4f vp = viewport();
  mat4f w2i = vp * p * v;    

  return w2i.inverse();
}


inline 
void Camera::frame(const float z, vec3f &dO, vec3f &dx, vec3f &dy) const
{
  mat4f i2w = image_toWorld(); // passage de l'image vers le monde

  vec3f uo(0, 0, z);
  vec3f ux(1, 0, z);
  vec3f uy(0, 1, z);
  dO = mult_vec_mat(uo, i2w);
  vec3f d1 = mult_vec_mat(ux, i2w); // récupère la colonne x
  vec3f d2 = mult_vec_mat(uy, i2w); // récupère la colonne y

  dx = subVec(dO, d1);
  dy = subVec(dO, d2);
}


struct NearestAxisResult
{
  vec3f nearestForwardAxis;
  vec3f nearestUpAxis;
};

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

  return {nearestForwardAxis, nearestUpAxis};
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
  vec3f forwardDirection = forward_direction();
  vec3f upDirection = up_direction();

  auto [nearestForwardAxis, nearestUpAxis] = nearest_axis(forwardDirection, upDirection);

  m_pitch = m_targetPitch;
  m_yaw = m_targetYaw;
  m_position = m_targetPosition;

  m_orientation = compute_rotation(nearestForwardAxis, nearestUpAxis);
}

#endif // CGAL_CAMERA_H
