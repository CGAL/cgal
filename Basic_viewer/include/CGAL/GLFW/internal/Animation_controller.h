#ifndef CGAL_ANIMATION_MANAGER_H
#define CGAL_ANIMATION_MANAGER_H

#include <vector>
#include <chrono>

#include "utils.h"

using namespace std::chrono_literals;

struct CameraKeyFrame 
{
  vec3f position;
  quatf orientation;

  CameraKeyFrame(const vec3f& _position, const quatf& _orientation) : 
    position(_position), 
    orientation(_orientation) 
  {
  }

  CameraKeyFrame() : 
    position(vec3f::Identity()), 
    orientation(quatf::Identity()) 
  {
  }
};

class Animation_controller 
{
public:
  using DurationType = std::chrono::milliseconds;
  using TimePoint = std::chrono::_V2::steady_clock::time_point;
  using KeyFrameBuffer = std::vector<CameraKeyFrame>;
  
public:
  void start();
  CameraKeyFrame run();
  void stop(const float frameNumber);

  void add_key_frame(const vec3f& position, const quatf& orientation);

  void set_duration(DurationType duration);

  inline bool is_running() const { return m_isRunning; }

  inline void clear_buffer() { m_keyFrames.clear(); }

  inline quatf get_rotation() const { return m_interpolatedRotation; }

  inline vec3f get_translation() const { return m_interpolatedTranslation; }

  inline float get_frame() const { return m_currentFrameNumber; }
  
private:
  CameraKeyFrame key_frame_interpolation(const float time);

  void compute_timestamp();

private: 
  quatf m_interpolatedRotation { quatf::Identity() };
  vec3f m_interpolatedTranslation { vec3f::Zero() };

  DurationType m_duration { 15s };   

  KeyFrameBuffer m_keyFrames {};

  TimePoint m_startTime {}; 

  CameraKeyFrame m_lastFrameData {};

  float m_lastFrameNumber    { 0.0f };   
  float m_currentFrameNumber { 0.0f };  
  float m_timestamp          { 0.0f };  

  bool m_isRunning { false };
};

/********************METHOD DEFINITIONS********************/

inline 
void Animation_controller::start() 
{
  if (m_keyFrames.size() <= 1) 
  {
    return;
  }

  if (!m_isRunning) 
  {
    m_startTime = std::chrono::steady_clock::now();
    m_isRunning = true;
  }
}

inline 
CameraKeyFrame Animation_controller::run() 
{
  auto now = std::chrono::steady_clock::now();
  float elapsedTime = static_cast<float>(std::chrono::duration_cast<std::chrono::milliseconds>(now - m_startTime).count());
  m_currentFrameNumber = m_lastFrameNumber + elapsedTime;
  if (m_currentFrameNumber >= static_cast<float>(m_duration.count())) 
  {
    m_currentFrameNumber = 0.0;
    stop(0.0);
    return m_lastFrameData;
  }

  m_lastFrameData = key_frame_interpolation(m_currentFrameNumber); 
  return m_lastFrameData;
}

inline 
void Animation_controller::stop(const float frameNumber) 
{
  m_isRunning = false;
  m_lastFrameNumber = frameNumber;
}

inline 
void Animation_controller::add_key_frame(const vec3f& position, const quatf& orientation) 
{ 
  CameraKeyFrame keyFrame(position, orientation);
  m_keyFrames.push_back(keyFrame);
  if (m_keyFrames.size() > 1) 
  {
    compute_timestamp();
  } 
}

inline 
void Animation_controller::set_duration(DurationType duration) 
{ 
  m_duration = duration; 
  if (m_keyFrames.size() > 1) 
  {
    compute_timestamp();
  }
}

static const float EPSILON = 0.001;

inline 
CameraKeyFrame Animation_controller::key_frame_interpolation(const float time)
{
  assert(m_timestamp > 0.0 + EPSILON || m_timestamp < 0.0 - EPSILON);

  float t = time / m_timestamp;

  unsigned int lowerIndex = std::floor(t);
  unsigned int upperIndex = std::ceil(t);

  CameraKeyFrame keyFrame0 = m_keyFrames.at(lowerIndex); 
  CameraKeyFrame keyFrame1 = m_keyFrames.at(upperIndex);

  t -= static_cast<float>(lowerIndex); 

  m_interpolatedRotation = keyFrame0.orientation.slerp(t, keyFrame1.orientation); 

  m_interpolatedTranslation = lerp(keyFrame0.position, keyFrame1.position, t);

  return {m_interpolatedTranslation, m_interpolatedRotation};
}

inline 
void Animation_controller::compute_timestamp() 
{ 
  assert(m_keyFrames.size() > 1);
  m_timestamp = static_cast<float>(m_duration.count()) / static_cast<float>(m_keyFrames.size() - 1); 
}

#endif // CGAL_ANIMATION_MANAGER_H
