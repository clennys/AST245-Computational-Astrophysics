#ifndef TIMER_H_
#define TIMER_H_

#include <chrono>
#include <vector>

using hr_clock = std::chrono::high_resolution_clock;

class Timer {
  public:
    std::chrono::time_point<hr_clock> m_start_point;
    std::chrono::time_point<hr_clock> m_end_point;
    std::vector<std::chrono::time_point<hr_clock>> m_checkpoints;
    long m_duration;

    Timer();
    auto start() -> void;
    auto stop() -> void;
    auto checkpoint() -> void;
    auto time_between_checkpoints(ulong i, ulong j) -> long;
    auto print_duration() -> void;
};

#endif // ! TIMER_H_
