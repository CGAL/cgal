
#ifndef MESSAGE_MANAGER_H
#define MESSAGE_MANAGER_H

#include <functional>
#include <map>
#include <string>


class Message_manager
{
public:
  static void add(const std::string& msg_name, std::function<void()> callback);
  static void notify_all(const std::string& msg_name);

protected:
  using Callbacks = std::vector<std::function<void()>>;
  static std::map<std::string, Callbacks>  s_message_map;

};


#endif
