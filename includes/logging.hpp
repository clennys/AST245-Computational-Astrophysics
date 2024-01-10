#ifndef LOG_H_
#define LOG_H_

#include <format>
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
static constexpr std::string_view RESET = "\033[0m";

} // namespace Colors

template <class... Args> auto dbg(std::format_string<Args...> fmt, Args &&...args) {
    std::cout << "[ DEBUG ]: " << std::format(fmt, std::forward<Args>(args)...) << "\n";
}

template <class... Args> auto info(std::format_string<Args...> fmt, Args &&...args) {
    std::cout << Colors::BRIGHT_BLUE << "[ INFO ]: " << Colors::RESET
              << std::format(fmt, std::forward<Args>(args)...) << "\n";
}

template <class... Args> auto warn(std::format_string<Args...> fmt, Args &&...args) {
    std::cerr << Colors::YELLOW << "[ WARN ]: " << Colors::RESET
              << std::format(fmt, std::forward<Args>(args)...) << "\n";
}

template <class... Args> auto err(std::format_string<Args...> fmt, Args &&...args) {
    std::cerr << Colors::RED << "[ ERROR ]: " << Colors::RESET
              << std::format(fmt, std::forward<Args>(args)...) << "\n";
}

} // namespace Logging

#endif // ! LOG_H_
