#include "timer.hpp"
#include <iostream>

Timer::Timer() {}

auto Timer::start() -> void {
    m_start_point = hr_clock::now();
    m_checkpoints.push_back(m_start_point);
}

auto Timer::stop() -> void {
    m_end_point = hr_clock::now();
    m_checkpoints.push_back(m_end_point);

    m_duration = time_between_checkpoints(0, m_checkpoints.size() - 1);
}

auto Timer::checkpoint() -> void { m_checkpoints.push_back(hr_clock::now()); }

auto Timer::print_duration() -> void {
    std::cout << "Time ellapsed: " << m_duration << "ms"
              << " -- " << static_cast<float>(m_duration) / 1000 << "s"
              << " -- " << static_cast<float>(m_duration) / 60000 << "m" << std::endl;
}

auto Timer::time_between_checkpoints(ulong i, ulong j) -> long {
    auto start = std::chrono::time_point_cast<std::chrono::milliseconds>(m_checkpoints[i])
                     .time_since_epoch()
                     .count();
    auto end = std::chrono::time_point_cast<std::chrono::milliseconds>(m_checkpoints[j])
                   .time_since_epoch()
                   .count();

    return end - start;
}
