#include "data.hpp"
#include "logging.hpp"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

auto Data::read_data(const std::filesystem::path &path_name)
    -> std::optional<std::vector<Particle3D>> {
    std::ifstream file(path_name);

    if (!file.is_open()) {
        Logging::err("Unable to open file");
        return std::nullopt;
    }

    std::string line;
    std::vector<Particle3D> particles;
    while (getline(file, line)) {
        std::stringstream ss(line);
        std::string word;
        // we have 10 values
        // Arrays no, Masses[i], x[i], y[i], z[i], Vx[i], Vy[i], Vz[i],
        // softening[i], potential[i]
        Particle3D part;
        int i = 0;

        while (ss >> word) {
            if (i == 0) {
                // we formally skip the array number
            } else if (i == 1) {
                part.mass = std::stod(word);
            } else if (i == 2) {
                part.position.x() = std::stod(word);
            } else if (i == 3) {
                part.position.y() = std::stod(word);
            } else if (i == 4) {
                part.position.z() = std::stod(word);
            } else if (i == 5) {
                part.velocity.x() = std::stod(word);
            } else if (i == 6) {
                part.velocity.y() = std::stod(word);
            } else if (i == 7) {
                part.velocity.z() = std::stod(word);
            } else if (i == 8) {
                part.softening = std::stod(word);
            } else if (i == 9) {
                part.potential = std::stod(word);
            }
            part.update_origin_dist();
            particles.push_back(part);
            i++;
            // part.print_summary();
        }
        // in case the file is formatted in an invalid way
        if (i != 10) {
            return std::nullopt;
        }
    }
    file.close();
    return particles;
}
