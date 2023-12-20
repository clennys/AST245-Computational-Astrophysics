#ifndef LOG_H_
#define LOG_H_

#include <iostream>
#include <string_view>

/**
 * TODO: (aver) add timing to loggin
 */

namespace Logging {
namespace Colors {

// Define some ANSI color codes
static constexpr std::string_view RED = "\033[31m";
static constexpr std::string_view YELLOW = "\033[33m";
static constexpr std::string_view BRIGHT_BLUE = "\033[94m";
static constexpr std::string_view RESET = "\033[0m";

} // namespace Colors

inline auto dbg(const std::string_view &dbg_msg) { std::cerr << "DEBUG: " << dbg_msg << "\n"; }

inline auto info(const std::string_view &info_msg) {
    std::cerr << Colors::BRIGHT_BLUE << "INFO: " << Colors::RESET << info_msg << "\n";
}

inline auto warn(const std::string_view &warn_msg) {
    std::cerr << Colors::YELLOW << "WARN: " << Colors::RESET << warn_msg << "\n";
}

inline auto err(const std::string_view &error_msg) {
    std::cerr << Colors::RED << "ERROR: " << Colors::RESET << error_msg << "\n";
}

} // namespace Logging

#endif // ! LOG_H_
