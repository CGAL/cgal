// Copyright (c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#include "Message_manager.h"

std::map<std::string, Message_manager::Callbacks>
Message_manager::s_message_map;

void Message_manager::add(const std::string& msg_name,
                          std::function<void()> callback)
{ s_message_map[msg_name].push_back(callback); }

void Message_manager::notify_all(const std::string& msg_name) {
  auto it = s_message_map.find(msg_name);
  if (s_message_map.cend() != it) {
    auto& callbacks = it->second;
    for (auto& cb : callbacks) cb();
  }
}
