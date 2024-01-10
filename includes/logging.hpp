#ifndef LOG_H_
#define LOG_H_

#include <chrono>
#include <format>
#include <iomanip>
#include <iostream>
#include <string_view>

/**
 * TODO: (aver)
 * - add timing to loggin
 */

namespace Logging {
namespace Colors {

// Define some ANSI color codes
static constexpr std::string_view RED = "\033[31m";
static constexpr std::string_view YELLOW = "\033[33m";
static constexpr std::string_view BRIGHT_BLUE = "\033[94m";
static constexpr std::string_view WHITE = "\033[37m";
static constexpr std::string_view RESET = "\033[0m";

} // namespace Colors

template <class... Args> auto dbg(std::format_string<Args...> fmt, Args &&...args) {
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    auto time = std::put_time(std::localtime(&in_time_t), "%X");

    std::cout << Colors::WHITE << "[ " << time << " DEBUG ]: " << Colors::RESET
              << std::format(fmt, std::forward<Args>(args)...) << "\n";
}

template <class... Args> auto info(std::format_string<Args...> fmt, Args &&...args) {
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    auto time = std::put_time(std::localtime(&in_time_t), "%X");

    std::cout << Colors::BRIGHT_BLUE << "[ " << time << "  INFO ]: " << Colors::RESET
              << std::format(fmt, std::forward<Args>(args)...) << "\n";
}

template <class... Args> auto warn(std::format_string<Args...> fmt, Args &&...args) {
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    auto time = std::put_time(std::localtime(&in_time_t), "%X");

    std::cerr << Colors::YELLOW << "[ " << time << "  WARN ]: " << Colors::RESET
              << std::format(fmt, std::forward<Args>(args)...) << "\n";
}

template <class... Args> auto err(std::format_string<Args...> fmt, Args &&...args) {
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    auto time = std::put_time(std::localtime(&in_time_t), "%X");

    std::cerr << Colors::RED << "[ " << time << " ERROR ]: " << Colors::RESET
              << std::format(fmt, std::forward<Args>(args)...) << "\n";
}

} // namespace Logging

#endif // ! LOG_H_
