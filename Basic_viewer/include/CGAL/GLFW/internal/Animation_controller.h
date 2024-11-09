#ifndef CGAL_GLFW_INTERNAL_ANIMATION_CONTROLLER_H
#define CGAL_GLFW_INTERNAL_ANIMATION_CONTROLLER_H

#include <vector>
#include <chrono>

#include "utils.h"

using namespace std::chrono_literals;

struct AnimationKeyFrame 
{
  vec3f position;
  quatf orientation;

  AnimationKeyFrame(const vec3f& _position, const quatf& _orientation) : 
    position(_position), 
    orientation(_orientation) 
  {
  }

  AnimationKeyFrame() : 
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
  using KeyFrameBuffer = std::vector<AnimationKeyFrame>;
  
public:
  void start();
  AnimationKeyFrame run();
  void stop(const float frameNumber);

  void add_key_frame(const vec3f& position, const quatf& orientation);

  void set_duration(DurationType duration);

  inline bool is_running() const { return m_IsRunning; }

  inline void clear_buffer() { m_KeyFrames.clear(); }

  inline quatf get_rotation() const { return m_InterpolatedRotation; }

  inline vec3f get_translation() const { return m_InterpolatedTranslation; }

  inline float get_frame() const { return m_CurrentFrameNumber; }

  inline size_t number_of_key_frames() const { return m_KeyFrames.size(); }
  
private:
  AnimationKeyFrame key_frame_interpolation(const float time);

  void compute_timestamp();

private: 
  quatf m_InterpolatedRotation { quatf::Identity() };
  vec3f m_InterpolatedTranslation { vec3f::Zero() };

  DurationType m_Duration { 5s };   

  KeyFrameBuffer m_KeyFrames {};

  TimePoint m_StartTime {}; 

  AnimationKeyFrame m_LastFrameData {};

  float m_LastFrameNumber    { 0.0f };   
  float m_CurrentFrameNumber { 0.0f };  
  float m_Timestamp          { 0.0f };  

  bool m_IsRunning { false };
};

/********************METHOD IMPLEMENTATIONS********************/

void Animation_controller::start() 
{
  if (m_KeyFrames.size() <= 1) 
  {
    return;
  }

  if (!m_IsRunning) 
  {
    m_StartTime = std::chrono::steady_clock::now();
    m_IsRunning = true;
  }
}

AnimationKeyFrame Animation_controller::run() 
{
  auto now = std::chrono::steady_clock::now();
  float elapsedTime = static_cast<float>(std::chrono::duration_cast<std::chrono::milliseconds>(now - m_StartTime).count());
  m_CurrentFrameNumber = m_LastFrameNumber + elapsedTime;
  if (m_CurrentFrameNumber >= static_cast<float>(m_Duration.count())) 
  {
    m_CurrentFrameNumber = 0.0;
    stop(0.0);
    return m_LastFrameData;
  }

  m_LastFrameData = key_frame_interpolation(m_CurrentFrameNumber); 
  return m_LastFrameData;
}

void Animation_controller::stop(const float frameNumber) 
{
  m_IsRunning = false;
  m_LastFrameNumber = frameNumber;
}

void Animation_controller::add_key_frame(const vec3f& position, const quatf& orientation) 
{ 
  AnimationKeyFrame keyFrame(position, orientation);
  m_KeyFrames.push_back(keyFrame);
  if (m_KeyFrames.size() > 1) 
  {
    compute_timestamp();
  } 
}

void Animation_controller::set_duration(DurationType duration) 
{ 
  m_Duration = duration; 
  if (m_KeyFrames.size() > 1) 
  {
    compute_timestamp();
  }
}

AnimationKeyFrame Animation_controller::key_frame_interpolation(const float time)
{
  assert(!utils::equal_float(m_Timestamp, 0));

  float t = time / m_Timestamp;

  unsigned int lowerIndex = std::floor(t);
  unsigned int upperIndex = std::ceil(t);

  AnimationKeyFrame keyFrame0 = m_KeyFrames.at(lowerIndex); 
  AnimationKeyFrame keyFrame1 = m_KeyFrames.at(upperIndex);

  t -= static_cast<float>(lowerIndex); 

  m_InterpolatedRotation = keyFrame0.orientation.slerp(t, keyFrame1.orientation); 

  m_InterpolatedTranslation = utils::lerp(keyFrame0.position, keyFrame1.position, t);

  return {m_InterpolatedTranslation, m_InterpolatedRotation};
}

void Animation_controller::compute_timestamp() 
{ 
  assert(m_KeyFrames.size() > 1);
  m_Timestamp = static_cast<float>(m_Duration.count()) / static_cast<float>(m_KeyFrames.size() - 1); 
}

#endif // CGAL_GLFW_INTERNAL_ANIMATION_CONTROLLER_H
