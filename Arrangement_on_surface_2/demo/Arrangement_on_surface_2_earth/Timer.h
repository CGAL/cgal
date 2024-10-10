// Copyright(c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <iostream>
#include <string>

class Timer {
public:
  Timer() { reset(); }

  void reset() { m_Start = std::chrono::high_resolution_clock::now(); }

  double elapsed() {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(
      std::chrono::high_resolution_clock::now() - m_Start).count() *
      0.001 * 0.001 * 0.001;
  }

  double elapsed_millis() { return elapsed() * 1000.0f; }

private:
  std::chrono::time_point<std::chrono::high_resolution_clock> m_Start;
};

class ScopedTimer {
public:
  ScopedTimer(const std::string& name) : m_name(name) {}

  ~ScopedTimer() {
    double time = m_timer.elapsed_millis();
    std::cout << "[TIMER] " << m_name << " - " << time << "ms\n";
  }

private:
  std::string m_name;
  Timer m_timer;
};


#endif
