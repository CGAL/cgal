#pragma once
#include <chrono>
#include <iostream>
#include "colorized_text.h"

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
        Colorized_print("[Time log] "+ _title + " > " + subTitle + " = " + std::to_string(static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(now - lastSubStep).count()) / 1000.) + "s. " + subText, ConsoleTextColor::Blue);
        lastSubStep = now;
    }

    void log_total_time() {
        auto now = std::chrono::steady_clock::now();
        Colorized_print("[Time log] "+ _title + " > " + " = " + std::to_string(static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(now - begin).count()) / 1000.) + "s. ", ConsoleTextColor::BrightBlue);
    }
private:
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point lastSubStep;

    std::string _title;
};

}

