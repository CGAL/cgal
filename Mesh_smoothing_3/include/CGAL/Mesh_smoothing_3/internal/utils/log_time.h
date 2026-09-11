// Copyright (c) 2026  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : François Protais

#ifndef CGAL_MESH_SMOOTHING_3_INTERNAL_UTILS_LOG_TIME_H
#define CGAL_MESH_SMOOTHING_3_INTERNAL_UTILS_LOG_TIME_H

#include <CGAL/license/Mesh_smoothing_3.h>

#include <CGAL/IO/Color_ostream.h>

#include <chrono>
#include <iostream>

namespace CGAL {

namespace Mesh_smoothing_3_internal {

class Time_log {
public:
    Time_log(std::string const &title): _title(title) {
        restart();
    }

    void restart() {
        begin = std::chrono::steady_clock::now();
        lastSubStep = begin;
    }

    void log_sub_step(std::string const &subTitle, std::string const &subText = "") {
        auto now = std::chrono::steady_clock::now();
        CGAL::IO::Color_stream_guard blue(std::cout, CGAL::IO::Ansi_color::Blue);
        std::cout << "[Time log] " << _title << " > " << subTitle << " = " << static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(now - lastSubStep).count()) / 1000. << "s. " << subText << std::endl;
        lastSubStep = now;
    }

    void log_total_time() {
        auto now = std::chrono::steady_clock::now();
        CGAL::IO::Color_stream_guard bright_blue(std::cout, CGAL::IO::Ansi_color::BrightBlue);
        std::cout << "[Time log] " << _title << " > " << " = " << static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(now - begin).count()) / 1000. << "s. " << std::endl;
    }
private:
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point lastSubStep;

    std::string _title;
};

} } // end of CGAL::Mesh_smoothing_3_internal namespace

#endif // CGAL_MESH_SMOOTHING_3_INTERNAL_UTILS_LOG_TIME_H
