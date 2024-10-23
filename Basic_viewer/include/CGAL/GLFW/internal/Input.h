#ifndef CGAL_INPUT_GLFW_H
#define CGAL_INPUT_GLFW_H

#include <GLFW/glfw3.h>
#include <unordered_map>
#include <vector>
#include <set>
#include <utility>
#include <string>
#include <cctype>
#include <iomanip>
#include <chrono>
#include <cassert>
#include <iostream>

using namespace std::chrono_literals;

enum class InputDevice { MOUSE, KEYBOARD };
enum class InputMode { HOLD, RELEASE };
enum class KeyboardLayout { QWERTY, AZERTY };

struct InputBinding 
{
  std::vector<int> keys {};
  InputMode mode {};
  InputDevice device {};

  inline int priority() const { return static_cast<int>(keys.size()); }  
  inline int size() const { return static_cast<int>(keys.size()); }  
  inline int get(int idx) const 
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

struct ActionComparator 
{
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

  void print_help();

  std::string get_binding_text_from_action(ActionEnum action);
  inline void set_action_description(const DescriptorType& descriptor) { m_Descriptor = descriptor; }

  inline std::pair<float, float> get_mouse_position() const { return m_MousePosition; }
  inline std::pair<float, float> get_mouse_delta() const { return m_MouseDelta; }
  inline std::pair<float, float> get_roll_mouse_delta() const { return m_RollMouseDelta; }
  
  inline float get_scroll_yOffset() { return clear_yOffset(); }
  inline void set_keyboard_layout(KeyboardLayout layout) { m_KeyboardLayout = layout; }
  inline bool is_qwerty_layout() const { return m_KeyboardLayout == KeyboardLayout::QWERTY; }

  bool is_key_pressed(GLFWwindow* window, int key) const;
  inline bool has_active_actions() const { return m_ActivatedActions.size() > 0; }

  static int map_to_azerty(int key);
  static int map_to_qwerty(int key);
  static std::string get_string_from_keycode(int key);

private:
  void register_action_events(const InputBinding& keys, ActionEnum action);


  inline float clear_yOffset() { int tmp = m_YOffset; m_YOffset = 0; return tmp; }

protected:
  virtual void start_action(ActionEnum action, const float deltaTime=0.0) = 0;
  virtual void action_event(ActionEnum action, const float deltaTime=0.0) = 0;
  virtual void end_action(ActionEnum action, const float deltaTime=0.0) = 0;

  virtual void double_click_event(int key) = 0;
  virtual void scroll_event(const float deltaTime) = 0;

  void on_key_event(int key, int scancode, int action, int mods);
  void on_cursor_event(double xpos, double ypos, int windowWidth, int windowHeight);
  void on_mouse_button_event(int button, int action, int mods);
  void on_scroll_event(double xoffset, double yoffset);

  void handle_events(const float deltaTime=0.0);
  
private:
  std::unordered_map<int, bool> m_KeyPressed {};
  std::unordered_map<int, bool> m_KeyHold {};
  std::unordered_map<int, bool> m_KeyConsumed {};
  
  std::unordered_map<ActionEnum, bool> m_StartedActions {};
  std::unordered_map<ActionEnum, bool> m_ActivatedActions {};

  DescriptorType m_Descriptor; 

  std::set<Action, ActionComparator> m_ActionEvents {};
  std::unordered_map<ActionEnum, std::vector<InputBinding>> m_ActionToBindingMapping {};

  std::pair<float, float> m_MousePosition {};
  std::pair<float, float> m_MouseDelta {};
  std::pair<float, float> m_RollMouseDelta {};

  KeyboardLayout m_KeyboardLayout { KeyboardLayout::QWERTY };

  float m_YOffset { 0 };

  std::chrono::_V2::system_clock::time_point m_LastTimeClick { std::chrono::high_resolution_clock::now() };
  int m_LastButtonClicked { GLFW_KEY_UNKNOWN };
};

void Input::add_keyboard_action(std::initializer_list<int> keys, InputMode mode, ActionEnum action)
{
  register_action_events(InputBinding {
    .keys = keys,
    .mode = mode,
    .device = InputDevice::KEYBOARD
  }, action);
}

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
  m_ActionEvents.insert(Action {
    .binding = binding,
    .action = action 
  });

  m_ActionToBindingMapping[action].emplace_back(binding);
}

void Input::on_key_event(int key, int scancode, int action, int mods)
{
  if (!is_qwerty_layout()) key = map_to_azerty(key);

  if (action == GLFW_PRESS || action == GLFW_REPEAT)
  {
    m_KeyPressed[key] = true;
    m_KeyHold[key] = true;
  }

  if (action == GLFW_RELEASE)
  {
    m_KeyHold.erase(key);
  }
}

void Input::on_cursor_event(double xpos, double ypos, int windowWidth, int windowHeight)
{
  auto& [mouseX, mouseY] = m_MousePosition; 
  m_MouseDelta = { xpos - mouseX, ypos - mouseY };
  m_MousePosition = { xpos, ypos };

  m_RollMouseDelta = m_MouseDelta;
  if (ypos < windowHeight / 2)
  {
    m_RollMouseDelta.first *= -1;
  }

  if (xpos < windowWidth / 2)
  {
    m_RollMouseDelta.second *= -1;
  }
}

void Input::on_scroll_event(double xoffset, double yoffset)
{
  m_YOffset = yoffset;
}

void Input::on_mouse_button_event(int button, int action, int mods)
{
  if (action == GLFW_PRESS) 
  {
    auto time = std::chrono::high_resolution_clock::now();

    const std::chrono::duration<double> dt = time - m_LastTimeClick;

    if (m_LastButtonClicked == button && dt < 250ms)
    {
      this->double_click_event(button);
    }

    m_LastTimeClick = time;
    m_LastButtonClicked = button;

    m_KeyPressed[button] = true;
    m_KeyHold[button] = true;
  }

  if (action == GLFW_RELEASE)
  {
    m_KeyHold.erase(button);
  }
}

void Input::handle_events(const float deltaTime)
{
  m_KeyPressed.clear();
  m_KeyConsumed.clear();
  m_ActivatedActions.clear();

  glfwPollEvents();

  if (m_YOffset != 0)
  {
    scroll_event(deltaTime);
  }

  for (auto& [binding, action] : m_ActionEvents)
  {
    std::unordered_map<int, bool> &mappedkeyState = binding.mode == InputMode::HOLD ? m_KeyHold : m_KeyPressed;

    if (!m_KeyConsumed[binding.get(0)] && mappedkeyState[binding.get(0)] && 
        (binding.priority() < 2 || m_KeyHold[binding.get(1)]) && 
        (binding.priority() < 3 || m_KeyHold[binding.get(2)])
      )
    {
      m_KeyConsumed[binding.get(0)] = true;

      m_ActivatedActions[action] = true;
      if (binding.mode == InputMode::HOLD)
      {
        if (!m_StartedActions[action])
        {
          m_StartedActions[action] = true;
          start_action(action, deltaTime);
        }
      }

      action_event(action, deltaTime);
    }
  }

  for (auto& [action, state] : m_StartedActions)
  {
    if (m_ActivatedActions.find(action) != m_ActivatedActions.end() && !m_ActivatedActions[action])
    {
      m_StartedActions[action] = false;
      end_action(action, deltaTime);
    }
  }

  m_MouseDelta = { 0.0f, 0.0f };
  m_RollMouseDelta = { 0.0f, 0.0f };
};

std::string Input::get_binding_text_from_action(ActionEnum action)  
{
  std::vector<InputBinding> bindings = m_ActionToBindingMapping[action];

  std::string result = "[";
  for (size_t i = 0; i < bindings.size(); ++i)
  {
    auto& binding = bindings[i];
    result += get_string_from_keycode(binding.get(binding.size()-1));
    for (int j = binding.size()-2; j >= 0; --j)
    {
      result = result + "+" + get_string_from_keycode(binding.get(j));
    }

    if (i < bindings.size()-1)
    {
      result += "] or [";
    }
  }

  result += "]";

  return result;
} 

void Input::add_description(ActionEnum action, const std::string& section, const std::string& description) 
{
  std::string binding = get_binding_text_from_action(action);
  add_description(binding, section, description);
}

void Input::add_description(const std::string& binding, const std::string& section, const std::string& description) 
{
  m_Descriptor[section].emplace_back(binding, description);
}

void Input::print_help()
{
  std::cout << "Basic Viewer GLFW - Features shortcuts  :" << "\n\n";

  int n = 140-2;
  std::string sectionFooter = "";
  for (int i = 0; i < n; ++i)
  {
    sectionFooter += "=";  
  }

  for (const auto& [sectionTitle, infos] : m_Descriptor) 
  {
    std::string sectionHeader = " Section : " + sectionTitle + " ";

    int nbIterations = n - sectionHeader.length();
    for (int i = 0; i < nbIterations/4; ++i)
    {
      sectionHeader = "=" + sectionHeader + "===";
    }
    std::cout << sectionHeader << "\n";

    for (const auto& [bindings, description] : infos) 
    {
      if (description.length() == 0) 
      {
        std::cout << '\n'; 
        continue;
      }

      std::cout << std::setw(n/4+4) << std::right  << bindings
                << " - "            << description << '\n';
    }
    
    std::cout << sectionFooter << "\n";
  }
  std::cout << std::endl;

  m_Descriptor.clear();
}

bool Input::is_key_pressed(GLFWwindow* window, int key) const
{
  if (m_KeyboardLayout == KeyboardLayout::AZERTY) 
  {
    key = map_to_qwerty(key);
  }

  int state = glfwGetKey(window, key);

  return state == GLFW_PRESS || state == GLFW_REPEAT;  
}

std::string Input::get_string_from_keycode(int key)
{
  if (key >= 39 && key <= 96) 
  {
    return std::string(1, static_cast<char>(key));
  }

  switch(key)
  {
    case GLFW_KEY_SPACE:           return "SPACE";
    case GLFW_KEY_ESCAPE:          return "ESCAPE";
    case GLFW_KEY_ENTER:           return "ENTER";
    case GLFW_KEY_TAB:             return "TAB";
    case GLFW_KEY_BACKSPACE:       return "BACKSPACE";
    case GLFW_KEY_INSERT:          return "INSERT";
    case GLFW_KEY_DELETE:          return "DELETE";
    case GLFW_KEY_RIGHT:           return "RIGHT";
    case GLFW_KEY_LEFT:            return "LEFT";
    case GLFW_KEY_DOWN:            return "DOWN";
    case GLFW_KEY_UP:              return "UP";
    case GLFW_KEY_PAGE_UP:         return "PAGE_UP";
    case GLFW_KEY_PAGE_DOWN:       return "PAGE_DOWN";
    case GLFW_KEY_HOME:            return "HOME";
    case GLFW_KEY_END:             return "END";
    case GLFW_KEY_CAPS_LOCK:       return "CAPS_LOCK";
    case GLFW_KEY_SCROLL_LOCK:     return "SCROLL_LOCK";
    case GLFW_KEY_NUM_LOCK:        return "NUM_LOCK";
    case GLFW_KEY_PRINT_SCREEN:    return "PRINT_SCREEN";
    case GLFW_KEY_PAUSE:           return "PAUSE";
    case GLFW_KEY_F1:              return "F1";
    case GLFW_KEY_F2:              return "F2";
    case GLFW_KEY_F3:              return "F3";
    case GLFW_KEY_F4:              return "F4";
    case GLFW_KEY_F5:              return "F5";
    case GLFW_KEY_F6:              return "F6";
    case GLFW_KEY_F7:              return "F7";
    case GLFW_KEY_F8:              return "F8";
    case GLFW_KEY_F9:              return "F9";
    case GLFW_KEY_F10:             return "F10";
    case GLFW_KEY_F11:             return "F11";
    case GLFW_KEY_F12:             return "F12";
    case GLFW_KEY_F13:             return "F13";
    case GLFW_KEY_F14:             return "F14";
    case GLFW_KEY_F15:             return "F15";
    case GLFW_KEY_F16:             return "F16";
    case GLFW_KEY_F17:             return "F17";
    case GLFW_KEY_F18:             return "F18";
    case GLFW_KEY_F19:             return "F19";
    case GLFW_KEY_F20:             return "F20";
    case GLFW_KEY_F21:             return "F21";
    case GLFW_KEY_F22:             return "F22";
    case GLFW_KEY_F23:             return "F23";
    case GLFW_KEY_F24:             return "F24";
    case GLFW_KEY_F25:             return "F25";
    case GLFW_KEY_KP_0:            return "NUMP_0";
    case GLFW_KEY_KP_1:            return "NUMP_1";
    case GLFW_KEY_KP_2:            return "NUMP_2";
    case GLFW_KEY_KP_3:            return "NUMP_3";
    case GLFW_KEY_KP_4:            return "NUMP_4";
    case GLFW_KEY_KP_5:            return "NUMP_5";
    case GLFW_KEY_KP_6:            return "NUMP_6";
    case GLFW_KEY_KP_7:            return "NUMP_7";
    case GLFW_KEY_KP_8:            return "NUMP_8";
    case GLFW_KEY_KP_9:            return "NUMP_9";
    case GLFW_KEY_KP_DECIMAL:      return "NUMP_DECIMAL";
    case GLFW_KEY_KP_DIVIDE:       return "NUMP_DIVIDE";
    case GLFW_KEY_KP_MULTIPLY:     return "NUMP_MULTIPLY";
    case GLFW_KEY_KP_SUBTRACT:     return "NUMP_SUBTRACT";
    case GLFW_KEY_KP_ADD:          return "NUMP_ADD";
    case GLFW_KEY_KP_ENTER:        return "NUMP_ENTER";
    case GLFW_KEY_KP_EQUAL:        return "NUMP_EQUAL";
    case GLFW_KEY_LEFT_SHIFT:      return "LSHIFT";
    case GLFW_KEY_LEFT_CONTROL:    return "LCTRL";
    case GLFW_KEY_LEFT_ALT:        return "LALT";
    case GLFW_KEY_LEFT_SUPER:      return "LSUPER";
    case GLFW_KEY_RIGHT_SHIFT:     return "RSHIFT";
    case GLFW_KEY_RIGHT_CONTROL:   return "RCTRL";
    case GLFW_KEY_RIGHT_ALT:       return "RALT";
    case GLFW_KEY_RIGHT_SUPER:     return "RSUPER";
    case GLFW_KEY_MENU:            return "MENU";  
    case GLFW_MOUSE_BUTTON_LAST:   return "MOUSE_BTN_LAST";
    case GLFW_MOUSE_BUTTON_LEFT:   return "MOUSE_BTN_LEFT";
    case GLFW_MOUSE_BUTTON_RIGHT:  return "MOUSE_BTN_RIGHT";
    case GLFW_MOUSE_BUTTON_MIDDLE: return "MOUSE_BTN_MIDDLE";         
    default:                       return "NONE";    
  }
}

int Input::map_to_azerty(int key)
{
  switch(key)
  {
    case GLFW_KEY_Q:         return GLFW_KEY_A;
    case GLFW_KEY_A:         return GLFW_KEY_Q;
    case GLFW_KEY_W:         return GLFW_KEY_Z;
    case GLFW_KEY_Z:         return GLFW_KEY_W;
    case GLFW_KEY_SEMICOLON: return GLFW_KEY_M;
    case GLFW_KEY_M:         return GLFW_KEY_COMMA;
    default:                 return key; 
  }
} 

int Input::map_to_qwerty(int key)
{
  switch(key)
  {
    case GLFW_KEY_Q:     return GLFW_KEY_A;
    case GLFW_KEY_A:     return GLFW_KEY_Q;
    case GLFW_KEY_W:     return GLFW_KEY_Z;
    case GLFW_KEY_Z:     return GLFW_KEY_W;
    case GLFW_KEY_M:     return GLFW_KEY_SEMICOLON;
    case GLFW_KEY_COMMA: return GLFW_KEY_M;
    default:             return key; 
  }
}


#endif // CGAL_INPUT_GLFW_H
