#ifndef LOG_H_
#define LOG_H_

#include <iostream>
#include <string_view>

/**
 * TODO: (aver) 
 * - add timing to loggin
 * - add std::format support (turned out to be more difficult than I though)
 */

namespace Logging {
namespace Colors {

// Define some ANSI color codes
static constexpr std::string_view RED = "\033[31m";
static constexpr std::string_view YELLOW = "\033[33m";
static constexpr std::string_view BRIGHT_BLUE = "\033[94m";
static constexpr std::string_view RESET = "\033[0m";

} // namespace Colors

template <class T> auto dbg(const T &dbg_msg) { std::cout << "[ DEBUG ]: " << dbg_msg << "\n"; }

template <class T> auto info(const T &info_msg) {
    std::cout << Colors::BRIGHT_BLUE << "[ INFO ]: " << Colors::RESET << info_msg << "\n";
}

template <class T> auto warn(const T &warn_msg) {
    std::cerr << Colors::YELLOW << "[ WARN ]: " << Colors::RESET << warn_msg << "\n";
}

template <class T> auto err(const T &error_msg) {
    std::cerr << Colors::RED << "[ ERROR ]: " << Colors::RESET << error_msg << "\n";
}

} // namespace Logging

#endif // ! LOG_H_
