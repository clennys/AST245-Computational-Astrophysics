#pragma once

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <boost/multiprecision/cpp_int.hpp>

/// @brief multiprecision integer
using mpint = boost::multiprecision::cpp_int;
/// implement formatting for mpint
template <>
struct fmt::formatter<mpint> : fmt::ostream_formatter {};
