// Copyright (c) 2023, 2024  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#ifndef MESSAGE_MANAGER_H
#define MESSAGE_MANAGER_H

#include <functional>
#include <map>
#include <string>
#include <vector>

class Message_manager {
public:
  static void add(const std::string& msg_name, std::function<void()> callback);
  static void notify_all(const std::string& msg_name);

protected:
  using Callbacks = std::vector<std::function<void()>>;
  static std::map<std::string, Callbacks> s_message_map;
};


#endif
