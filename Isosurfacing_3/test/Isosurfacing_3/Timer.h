#ifndef CGAL_TIMER_H
#define CGAL_TIMER_H

#include <chrono>
#include <iostream>

class ScopeTimer {
    using Clock = std::chrono::steady_clock;
    using TimePoint = Clock::time_point;

public:
    ScopeTimer() : start(Clock::now()), msg("Duration"), running(true) {}

    explicit ScopeTimer(const std::string& msg) : start(Clock::now()), msg(msg), running(true) {
        std::cout << msg << "..." << std::endl;
    }

    int64_t stop() {
        TimePoint end = Clock::now();
        running = false;
        return std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    }

    ~ScopeTimer() {
        if (running) {
            TimePoint end = Clock::now();
            int64_t duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            std::cout << msg << ": " << duration << " ms" << std::endl;
        }
    }

private:
    TimePoint start;
    std::string msg;
    bool running;
};

#endif  // CGAL_TIMER_H
