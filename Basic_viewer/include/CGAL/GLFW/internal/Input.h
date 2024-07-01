#ifndef CGAL_INPUT_GLFW_H
#define CGAL_INPUT_GLFW_H

#include <GLFW/glfw3.h>
#include <unordered_map>
#include <vector>
#include <set>
#include <string>
#include <cctype>
#include <iomanip>
#include <eigen3/Eigen/Core>
#include <chrono>

#include "key_code.h"

using namespace std::chrono_literals;

/**
 * TODO: Peux mieux faire, utiliser (int mod) pour shift/alt
 */

enum class InputDevice { MOUSE, KEYBOARD };
enum class InputMode { HOLD, RELEASE };
enum class KeyboardLayout { QWERTY, AZERTY };

struct InputBinding 
{
  std::vector<int> keys {};
  InputMode mode {};
  InputDevice device {};

  int priority() const { return static_cast<int>(keys.size()); }  
  int get(int idx) const 
  { 
    assert(static_cast<size_t>(idx) < keys.size()); 
    return keys.at(idx); 
  } 
};

struct Action
{
  InputBinding binding {};
  int action {};
};

struct ActionComparator {
  bool operator()(const Action& action_1, const Action& action_2) const 
  {

    return action_1.binding.priority() >= action_2.binding.priority();
  }
};

class Input
{
public:
  using ActionEnum = int;
  using DescriptorType = std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>>;
public:
  void add_keyboard_action(std::initializer_list<int> keys, InputMode mode, ActionEnum action); 

  void add_mouse_action(std::initializer_list<int> keys, InputMode mode, ActionEnum action); 

  void add_description(ActionEnum action, const std::string& section, const std::string& description);
  void add_description(const std::string& action, const std::string& section, const std::string& description);

  void print_help() const;

  std::string get_binding_text_from_action(ActionEnum action);
  inline void set_action_description(const DescriptorType& descriptor) { m_descriptor = descriptor; }

  inline Eigen::Vector2f get_mouse_position() const { return m_mousePosition; }
  inline Eigen::Vector2f get_mouse_delta() const { return m_mouseDelta; }
  
  inline float get_scroll_yOffset() { return clear_yOffset(); }
  inline void set_keyboard_layout(KeyboardLayout layout) { m_keyboardLayout = layout; }
  inline bool is_qwerty_layout() const { return m_keyboardLayout == KeyboardLayout::QWERTY; }

  static int map_to_azerty(int key);
  static int get_int_from_key(int key);
  static std::string get_string_from_keycode(int key);

private:
  void register_action_events(const InputBinding& keys, ActionEnum action);


  inline float clear_yOffset() { int tmp = m_yOffset; m_yOffset = 0; return tmp; }

protected:
  virtual void start_action(ActionEnum action, const double deltaTime=0.0) = 0;
  virtual void action_event(ActionEnum action, const double deltaTime=0.0) = 0;
  virtual void end_action(ActionEnum action, const double deltaTime=0.0) = 0;

  virtual void double_click_event(int key) = 0;
  virtual void scroll_event(const double deltaTime) = 0;

  void on_key_event(int key, int scancode, int action, int mods);
  void on_cursor_event(double xpos, double ypos);
  void on_mouse_button_event(int button, int action, int mods);
  void on_scroll_event(double xoffset, double yoffset);

  void handle_events(const double deltaTime=0.00);
  
private:
  std::unordered_map<int, bool> m_keyPressed {};
  std::unordered_map<int, bool> m_keyHold {};
  std::unordered_map<int, bool> m_keyConsumed {};
  
  std::unordered_map<ActionEnum, bool> m_startedActions {};
  std::unordered_map<ActionEnum, bool> m_activatedActions {};

  DescriptorType m_descriptor; 

  std::set<Action, ActionComparator> m_actionEvents {};
  std::unordered_map<ActionEnum, InputBinding> m_actionToBindingMapping {};

  Eigen::Vector2f m_mousePosition { Eigen::Vector2f::Zero() };
  Eigen::Vector2f m_mouseDelta { Eigen::Vector2f::Zero() };

  KeyboardLayout m_keyboardLayout { KeyboardLayout::QWERTY };

  float m_yOffset { 0 };

  std::chrono::_V2::system_clock::time_point m_lastClickTime { std::chrono::high_resolution_clock::now() };
  int m_lastButtonClicked { CGAL_BV_KEY_UNKNOWN };
};

inline 
void Input::add_keyboard_action(std::initializer_list<int> keys, InputMode mode, ActionEnum action)
{
  register_action_events(InputBinding {
    .keys = keys,
    .mode = mode,
    .device = InputDevice::KEYBOARD
  }, action);
}

inline 
void Input::add_mouse_action(std::initializer_list<int> keys, InputMode mode, ActionEnum action)
{
  register_action_events(InputBinding {
    .keys = keys,
    .mode = mode,
    .device = InputDevice::MOUSE
  }, action);
}

void Input::register_action_events(const InputBinding& binding, ActionEnum action)
{
  m_actionEvents.insert(Action {
    .binding = binding,
    .action = action 
  });

  m_actionToBindingMapping[action] = binding;
}

void Input::on_key_event(int key, int scancode, int action, int mods)
{
  if (!is_qwerty_layout()) key = map_to_azerty(key);
  int iKey = get_int_from_key(key);

  if (action == GLFW_PRESS || action == GLFW_REPEAT)
  {
    m_keyPressed[iKey] = true;
    m_keyHold[iKey] = true;
  }

  if (action == GLFW_RELEASE)
  {
    m_keyHold.erase(iKey);
  }
}

void Input::on_cursor_event(double xpos, double ypos)
{
  m_mouseDelta << 
    xpos - m_mousePosition.x(), 
    ypos - m_mousePosition.y();

  m_mousePosition << xpos, ypos;
}

void Input::on_scroll_event(double xoffset, double yoffset)
{
  m_yOffset = yoffset;
}

void Input::on_mouse_button_event(int button, int action, int mods)
{
  int iKey = get_int_from_key(button);

  if (action == GLFW_PRESS) 
  {
    auto time = std::chrono::high_resolution_clock::now();

    const std::chrono::duration<double> dt = time - m_lastClickTime;

    if (m_lastButtonClicked == iKey && dt < 250ms)
    {
      this->double_click_event(iKey);
    }

    m_lastClickTime = time;
    m_lastButtonClicked = iKey;

    m_keyPressed[iKey] = true;
    m_keyHold[iKey] = true;
  }

  if (action == GLFW_RELEASE)
  {
    m_keyHold.erase(iKey);
  }
}

void Input::handle_events(const double deltaTime)
{
  m_keyPressed.clear();
  m_keyConsumed.clear();
  m_activatedActions.clear();

  glfwPollEvents();

  if (m_yOffset != 0)
  {
    this->scroll_event(deltaTime);
  }

  for (auto& [binding, action] : m_actionEvents)
  {
    std::unordered_map<int, bool> &mappedkeyState = binding.mode == InputMode::HOLD ? m_keyHold : m_keyPressed;

    if (!m_keyConsumed[binding.get(0)] && mappedkeyState[binding.get(0)] && 
        (binding.priority() < 2 || m_keyHold[binding.get(1)]) && 
        (binding.priority() < 3 || m_keyHold[binding.get(2)])
      )
    {
      m_keyConsumed[binding.get(0)] = true;

      if (binding.mode == InputMode::HOLD)
      {
        m_activatedActions[action] = true;
        if (!m_startedActions[action])
        {
          m_startedActions[action] = true;
          start_action(action, deltaTime);
        }
      }

      this->action_event(action, deltaTime);
    }
  }

  for (auto& [action, state] : m_startedActions)
  {
    if (!m_activatedActions[action])
    {
      m_startedActions[action] = false;
      end_action(action, deltaTime);
    }
  }

  m_mouseDelta << 0.0f, 0.0f;
};

inline 
std::string Input::get_binding_text_from_action(ActionEnum action)  
{
  InputBinding binding = m_actionToBindingMapping[action];
  std::string ret = get_string_from_keycode(binding.get(0));
  for (size_t i = 1; i < binding.priority(); ++i)
  {
    ret = ret + "+" + get_string_from_keycode(binding.get(i));
  }

  return ret;
} 

inline
void Input::add_description(ActionEnum action, const std::string& section, const std::string& description) 
{
  std::string binding = get_binding_text_from_action(action);
  add_description(binding, section, description);
}

inline
void Input::add_description(const std::string& binding, const std::string& section, const std::string& description) 
{
  m_descriptor[section].emplace_back(binding, description);
}

void Input::print_help() const
{
  std::cout << "Basic Viewer GLFW - Features shortcute :" << "\n\n";

  int n = 100-2;
  std::string sectionFooter = "";
  for (int i = 0; i < n; ++i)
  {
    sectionFooter += "=";  
  }

  for (const auto& [sectionTitle, infos] : m_descriptor) 
  {
    std::string sectionHeader = " Section : " + sectionTitle + " ";

    int nbIterations = n - sectionHeader.length();
    for (int i = 0; i < nbIterations/4; ++i)
    {
      sectionHeader = "=" + sectionHeader + "===";
    }
    std::cout << sectionHeader << "\n";

    for (const auto& [binding, description] : infos) 
    {
      if (description.length() == 0) 
      {
        std::cout << '\n'; 
        continue;
      }

      std::cout << std::setw(n/4+4) << std::right << "[" + binding + "]"
                << " - " << description << '\n';
    }
    
    std::cout << sectionFooter << "\n";
  }
  std::cout << std::endl;
}

inline
std::string Input::get_string_from_keycode(int key)
{
  if (key >= 32 && key <= 96) 
  {
    return std::string(1, static_cast<char>(key));
  }

  switch(key)
  {
    case CGAL_BV_KEY_ESCAPE:                  return "ESCAPE";
    case CGAL_BV_KEY_ENTER:                   return "ENTER";
    case CGAL_BV_KEY_TAB:                     return "TAB";
    case CGAL_BV_KEY_BACKSPACE:               return "BACKSPACE";
    case CGAL_BV_KEY_INSERT:                  return "INSERT";
    case CGAL_BV_KEY_DELETE:                  return "DELETE";
    case CGAL_BV_KEY_RIGHT:                   return "RIGHT";
    case CGAL_BV_KEY_LEFT:                    return "LEFT";
    case CGAL_BV_KEY_DOWN:                    return "DOWN";
    case CGAL_BV_KEY_UP:                      return "UP";
    case CGAL_BV_KEY_PAGE_UP:                 return "PAGE_UP";
    case CGAL_BV_KEY_PAGE_DOWN:               return "PAGE_DOWN";
    case CGAL_BV_KEY_HOME:                    return "HOME";
    case CGAL_BV_KEY_END:                     return "END";
    case CGAL_BV_KEY_CAPS_LOCK:               return "CAPS_LOCK";
    case CGAL_BV_KEY_SCROLL_LOCK:             return "SCROLL_LOCK";
    case CGAL_BV_KEY_NUM_LOCK:                return "NUM_LOCK";
    case CGAL_BV_KEY_PRINT_SCREEN:            return "PRINT_SCREEN";
    case CGAL_BV_KEY_PAUSE:                   return "PAUSE";
    case CGAL_BV_KEY_F1:                      return "F1";
    case CGAL_BV_KEY_F2:                      return "F2";
    case CGAL_BV_KEY_F3:                      return "F3";
    case CGAL_BV_KEY_F4:                      return "F4";
    case CGAL_BV_KEY_F5:                      return "F5";
    case CGAL_BV_KEY_F6:                      return "F6";
    case CGAL_BV_KEY_F7:                      return "F7";
    case CGAL_BV_KEY_F8:                      return "F8";
    case CGAL_BV_KEY_F9:                      return "F9";
    case CGAL_BV_KEY_F10:                     return "F10";
    case CGAL_BV_KEY_F11:                     return "F11";
    case CGAL_BV_KEY_F12:                     return "F12";
    case CGAL_BV_KEY_F13:                     return "F13";
    case CGAL_BV_KEY_F14:                     return "F14";
    case CGAL_BV_KEY_F15:                     return "F15";
    case CGAL_BV_KEY_F16:                     return "F16";
    case CGAL_BV_KEY_F17:                     return "F17";
    case CGAL_BV_KEY_F18:                     return "F18";
    case CGAL_BV_KEY_F19:                     return "F19";
    case CGAL_BV_KEY_F20:                     return "F20";
    case CGAL_BV_KEY_F21:                     return "F21";
    case CGAL_BV_KEY_F22:                     return "F22";
    case CGAL_BV_KEY_F23:                     return "F23";
    case CGAL_BV_KEY_F24:                     return "F24";
    case CGAL_BV_KEY_F25:                     return "F25";
    case CGAL_BV_KEY_KP_0:                    return "NUMP_0";
    case CGAL_BV_KEY_KP_1:                    return "NUMP_1";
    case CGAL_BV_KEY_KP_2:                    return "NUMP_2";
    case CGAL_BV_KEY_KP_3:                    return "NUMP_3";
    case CGAL_BV_KEY_KP_4:                    return "NUMP_4";
    case CGAL_BV_KEY_KP_5:                    return "NUMP_5";
    case CGAL_BV_KEY_KP_6:                    return "NUMP_6";
    case CGAL_BV_KEY_KP_7:                    return "NUMP_7";
    case CGAL_BV_KEY_KP_8:                    return "NUMP_8";
    case CGAL_BV_KEY_KP_9:                    return "NUMP_9";
    case CGAL_BV_KEY_KP_DECIMAL:              return "NUMP_DECIMAL";
    case CGAL_BV_KEY_KP_DIVIDE:               return "NUMP_DIVIDE";
    case CGAL_BV_KEY_KP_MULTIPLY:             return "NUMP_MULTIPLY";
    case CGAL_BV_KEY_KP_SUBTRACT:             return "NUMP_SUBTRACT";
    case CGAL_BV_KEY_KP_ADD:                  return "NUMP_ADD";
    case CGAL_BV_KEY_KP_ENTER:                return "NUMP_ENTER";
    case CGAL_BV_KEY_KP_EQUAL:                return "NUMP_EQUAL";
    case CGAL_BV_KEY_LEFT_SHIFT:              return "LSHIFT";
    case CGAL_BV_KEY_LEFT_CONTROL:            return "LCTRL";
    case CGAL_BV_KEY_LEFT_ALT:                return "LALT";
    case CGAL_BV_KEY_LEFT_SUPER:              return "LSUPER";
    case CGAL_BV_KEY_RIGHT_SHIFT:             return "RSHIFT";
    case CGAL_BV_KEY_RIGHT_CONTROL:           return "RCTRL";
    case CGAL_BV_KEY_RIGHT_ALT:               return "RALT";
    case CGAL_BV_KEY_RIGHT_SUPER:             return "RSUPER";
    case CGAL_BV_KEY_MENU:                    return "MENU";  
    case CGAL_BV_MOUSE_BUTTON_LAST:           return "MOUSE_BTN_LAST";
    case CGAL_BV_MOUSE_BUTTON_LEFT:           return "MOUSE_BTN_LEFT";
    case CGAL_BV_MOUSE_BUTTON_RIGHT:          return "MOUSE_BTN_RIGHT";
    case CGAL_BV_MOUSE_BUTTON_MIDDLE:         return "MOUSE_BTN_MIDDLE";         
    default:                                  return "NONE";    
  }
}

inline 
int Input::map_to_azerty(int key)
{
  switch(key)
  {
    case GLFW_KEY_Q:                    return GLFW_KEY_A;
    case GLFW_KEY_A:                    return GLFW_KEY_Q;
    case GLFW_KEY_W:                    return GLFW_KEY_Z;
    case GLFW_KEY_Z:                    return GLFW_KEY_W;
    case GLFW_KEY_SEMICOLON:            return GLFW_KEY_M;
    case GLFW_KEY_M:                    return GLFW_KEY_COMMA;
    default:                            return key; 
  }
}


inline 
int Input::get_int_from_key(const int key)
{
  switch(key)
  {
    case GLFW_KEY_SPACE:                return CGAL_BV_KEY_SPACE;        
    case GLFW_KEY_APOSTROPHE:           return CGAL_BV_KEY_APOSTROPHE;   
    case GLFW_KEY_COMMA:                return CGAL_BV_KEY_COMMA;        
    case GLFW_KEY_MINUS:                return CGAL_BV_KEY_MINUS;        
    case GLFW_KEY_PERIOD:               return CGAL_BV_KEY_PERIOD;       
    case GLFW_KEY_SLASH:                return CGAL_BV_KEY_SLASH;        
    case GLFW_KEY_0:                    return CGAL_BV_KEY_0;            
    case GLFW_KEY_1:                    return CGAL_BV_KEY_1;            
    case GLFW_KEY_2:                    return CGAL_BV_KEY_2;            
    case GLFW_KEY_3:                    return CGAL_BV_KEY_3;            
    case GLFW_KEY_4:                    return CGAL_BV_KEY_4;            
    case GLFW_KEY_5:                    return CGAL_BV_KEY_5;            
    case GLFW_KEY_6:                    return CGAL_BV_KEY_6;            
    case GLFW_KEY_7:                    return CGAL_BV_KEY_7;            
    case GLFW_KEY_8:                    return CGAL_BV_KEY_8;            
    case GLFW_KEY_9:                    return CGAL_BV_KEY_9;            
    case GLFW_KEY_SEMICOLON:            return CGAL_BV_KEY_SEMICOLON;    
    case GLFW_KEY_EQUAL:                return CGAL_BV_KEY_EQUAL;        
    case GLFW_KEY_A:                    return CGAL_BV_KEY_A;            
    case GLFW_KEY_B:                    return CGAL_BV_KEY_B;            
    case GLFW_KEY_C:                    return CGAL_BV_KEY_C;            
    case GLFW_KEY_D:                    return CGAL_BV_KEY_D;            
    case GLFW_KEY_E:                    return CGAL_BV_KEY_E;            
    case GLFW_KEY_F:                    return CGAL_BV_KEY_F;            
    case GLFW_KEY_G:                    return CGAL_BV_KEY_G;            
    case GLFW_KEY_H:                    return CGAL_BV_KEY_H;            
    case GLFW_KEY_I:                    return CGAL_BV_KEY_I;            
    case GLFW_KEY_J:                    return CGAL_BV_KEY_J;            
    case GLFW_KEY_K:                    return CGAL_BV_KEY_K;            
    case GLFW_KEY_L:                    return CGAL_BV_KEY_L;            
    case GLFW_KEY_M:                    return CGAL_BV_KEY_M;            
    case GLFW_KEY_N:                    return CGAL_BV_KEY_N;            
    case GLFW_KEY_O:                    return CGAL_BV_KEY_O;            
    case GLFW_KEY_P:                    return CGAL_BV_KEY_P;            
    case GLFW_KEY_Q:                    return CGAL_BV_KEY_Q;            
    case GLFW_KEY_R:                    return CGAL_BV_KEY_R;            
    case GLFW_KEY_S:                    return CGAL_BV_KEY_S;            
    case GLFW_KEY_T:                    return CGAL_BV_KEY_T;            
    case GLFW_KEY_U:                    return CGAL_BV_KEY_U;            
    case GLFW_KEY_V:                    return CGAL_BV_KEY_V;            
    case GLFW_KEY_W:                    return CGAL_BV_KEY_W;            
    case GLFW_KEY_X:                    return CGAL_BV_KEY_X;            
    case GLFW_KEY_Y:                    return CGAL_BV_KEY_Y;            
    case GLFW_KEY_Z:                    return CGAL_BV_KEY_Z;            
    case GLFW_KEY_LEFT_BRACKET:         return CGAL_BV_KEY_LEFT_BRACKET; 
    case GLFW_KEY_BACKSLASH:            return CGAL_BV_KEY_BACKSLASH;    
    case GLFW_KEY_RIGHT_BRACKET:        return CGAL_BV_KEY_RIGHT_BRACKET;
    case GLFW_KEY_GRAVE_ACCENT:         return CGAL_BV_KEY_GRAVE_ACCENT; 
    case GLFW_KEY_WORLD_1:              return CGAL_BV_KEY_WORLD_1;      
    case GLFW_KEY_WORLD_2:              return CGAL_BV_KEY_WORLD_2;      
    
    case GLFW_KEY_ESCAPE:               return CGAL_BV_KEY_ESCAPE;                
    case GLFW_KEY_ENTER:                return CGAL_BV_KEY_ENTER;                
    case GLFW_KEY_TAB:                  return CGAL_BV_KEY_TAB;                  
    case GLFW_KEY_BACKSPACE:            return CGAL_BV_KEY_BACKSPACE;            
    case GLFW_KEY_INSERT:               return CGAL_BV_KEY_INSERT;                
    case GLFW_KEY_DELETE:               return CGAL_BV_KEY_DELETE;                
    case GLFW_KEY_RIGHT:                return CGAL_BV_KEY_RIGHT;                
    case GLFW_KEY_LEFT:                 return CGAL_BV_KEY_LEFT;                  
    case GLFW_KEY_DOWN:                 return CGAL_BV_KEY_DOWN;                  
    case GLFW_KEY_UP:                   return CGAL_BV_KEY_UP;                    
    case GLFW_KEY_PAGE_UP:              return CGAL_BV_KEY_PAGE_UP;              
    case GLFW_KEY_PAGE_DOWN:            return CGAL_BV_KEY_PAGE_DOWN;            
    case GLFW_KEY_HOME:                 return CGAL_BV_KEY_HOME;                  
    case GLFW_KEY_END:                  return CGAL_BV_KEY_END;                  
    case GLFW_KEY_CAPS_LOCK:            return CGAL_BV_KEY_CAPS_LOCK;            
    case GLFW_KEY_SCROLL_LOCK:          return CGAL_BV_KEY_SCROLL_LOCK;          
    case GLFW_KEY_NUM_LOCK:             return CGAL_BV_KEY_NUM_LOCK;              
    case GLFW_KEY_PRINT_SCREEN:         return CGAL_BV_KEY_PRINT_SCREEN;          
    case GLFW_KEY_PAUSE:                return CGAL_BV_KEY_PAUSE;                
    case GLFW_KEY_F1:                   return CGAL_BV_KEY_F1;                    
    case GLFW_KEY_F2:                   return CGAL_BV_KEY_F2;                    
    case GLFW_KEY_F3:                   return CGAL_BV_KEY_F3;                    
    case GLFW_KEY_F4:                   return CGAL_BV_KEY_F4;                    
    case GLFW_KEY_F5:                   return CGAL_BV_KEY_F5;                    
    case GLFW_KEY_F6:                   return CGAL_BV_KEY_F6;                    
    case GLFW_KEY_F7:                   return CGAL_BV_KEY_F7;                    
    case GLFW_KEY_F8:                   return CGAL_BV_KEY_F8;                    
    case GLFW_KEY_F9:                   return CGAL_BV_KEY_F9;                    
    case GLFW_KEY_F10:                  return CGAL_BV_KEY_F10;                  
    case GLFW_KEY_F11:                  return CGAL_BV_KEY_F11;                  
    case GLFW_KEY_F12:                  return CGAL_BV_KEY_F12;                  
    case GLFW_KEY_F13:                  return CGAL_BV_KEY_F13;                  
    case GLFW_KEY_F14:                  return CGAL_BV_KEY_F14;                  
    case GLFW_KEY_F15:                  return CGAL_BV_KEY_F15;                  
    case GLFW_KEY_F16:                  return CGAL_BV_KEY_F16;                  
    case GLFW_KEY_F17:                  return CGAL_BV_KEY_F17;                  
    case GLFW_KEY_F18:                  return CGAL_BV_KEY_F18;                  
    case GLFW_KEY_F19:                  return CGAL_BV_KEY_F19;                  
    case GLFW_KEY_F20:                  return CGAL_BV_KEY_F20;                  
    case GLFW_KEY_F21:                  return CGAL_BV_KEY_F21;                  
    case GLFW_KEY_F22:                  return CGAL_BV_KEY_F22;                  
    case GLFW_KEY_F23:                  return CGAL_BV_KEY_F23;                  
    case GLFW_KEY_F24:                  return CGAL_BV_KEY_F24;                  
    case GLFW_KEY_F25:                  return CGAL_BV_KEY_F25;                  
    case GLFW_KEY_KP_0:                 return CGAL_BV_KEY_KP_0;                  
    case GLFW_KEY_KP_1:                 return CGAL_BV_KEY_KP_1;                  
    case GLFW_KEY_KP_2:                 return CGAL_BV_KEY_KP_2;                  
    case GLFW_KEY_KP_3:                 return CGAL_BV_KEY_KP_3;                  
    case GLFW_KEY_KP_4:                 return CGAL_BV_KEY_KP_4;                  
    case GLFW_KEY_KP_5:                 return CGAL_BV_KEY_KP_5;                  
    case GLFW_KEY_KP_6:                 return CGAL_BV_KEY_KP_6;                  
    case GLFW_KEY_KP_7:                 return CGAL_BV_KEY_KP_7;                  
    case GLFW_KEY_KP_8:                 return CGAL_BV_KEY_KP_8;                  
    case GLFW_KEY_KP_9:                 return CGAL_BV_KEY_KP_9;                  
    case GLFW_KEY_KP_DECIMAL:           return CGAL_BV_KEY_KP_DECIMAL;            
    case GLFW_KEY_KP_DIVIDE:            return CGAL_BV_KEY_KP_DIVIDE;            
    case GLFW_KEY_KP_MULTIPLY:          return CGAL_BV_KEY_KP_MULTIPLY;          
    case GLFW_KEY_KP_SUBTRACT:          return CGAL_BV_KEY_KP_SUBTRACT;          
    case GLFW_KEY_KP_ADD:               return CGAL_BV_KEY_KP_ADD;                
    case GLFW_KEY_KP_ENTER:             return CGAL_BV_KEY_KP_ENTER;              
    case GLFW_KEY_KP_EQUAL:             return CGAL_BV_KEY_KP_EQUAL;              
    case GLFW_KEY_LEFT_SHIFT:           return CGAL_BV_KEY_LEFT_SHIFT;            
    case GLFW_KEY_LEFT_CONTROL:         return CGAL_BV_KEY_LEFT_CONTROL;          
    case GLFW_KEY_LEFT_ALT:             return CGAL_BV_KEY_LEFT_ALT;              
    case GLFW_KEY_LEFT_SUPER:           return CGAL_BV_KEY_LEFT_SUPER;            
    case GLFW_KEY_RIGHT_SHIFT:          return CGAL_BV_KEY_RIGHT_SHIFT;          
    case GLFW_KEY_RIGHT_CONTROL:        return CGAL_BV_KEY_RIGHT_CONTROL;        
    case GLFW_KEY_RIGHT_ALT:            return CGAL_BV_KEY_RIGHT_ALT;            
    case GLFW_KEY_RIGHT_SUPER:          return CGAL_BV_KEY_RIGHT_SUPER;          
    case GLFW_KEY_MENU:                 return CGAL_BV_KEY_MENU;                  

    case GLFW_MOUSE_BUTTON_LAST:        return CGAL_BV_MOUSE_BUTTON_LAST;    
    case GLFW_MOUSE_BUTTON_LEFT:        return CGAL_BV_MOUSE_BUTTON_LEFT;    
    case GLFW_MOUSE_BUTTON_RIGHT:       return CGAL_BV_MOUSE_BUTTON_RIGHT;    
    case GLFW_MOUSE_BUTTON_MIDDLE:      return CGAL_BV_MOUSE_BUTTON_MIDDLE;  

    default:                            return CGAL_BV_KEY_UNKNOWN;
  }
} 

#endif // CGAL_INPUT_GLFW_H
