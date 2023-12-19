#ifndef DATA_H_
#define DATA_H_

#include <filesystem>
#include <optional>
#include <vector>

#include "particle.hpp"

namespace Data {

/// @brief Reads in data for particles from a file path and returns a parsed list, null-opt if
/// failed
auto read_data(const std::filesystem::path &path_name) -> std::optional<std::vector<Particle3D>>;

} // namespace Data

#endif // ! DATA_H_
