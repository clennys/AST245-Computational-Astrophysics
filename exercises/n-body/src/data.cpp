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
                part.m_id = std::stoul(word);
            } else if (i == 1) {
                part.m_mass = std::stod(word);
            } else if (i == 2) {
                part.m_position.x() = std::stod(word);
            } else if (i == 3) {
                part.m_position.y() = std::stod(word);
            } else if (i == 4) {
                part.m_position.z() = std::stod(word);
            } else if (i == 5) {
                part.m_velocity.x() = std::stod(word);
            } else if (i == 6) {
                part.m_velocity.y() = std::stod(word);
            } else if (i == 7) {
                part.m_velocity.z() = std::stod(word);
            } else if (i == 8) {
                // WARN: (aver) ignore softening for now
                // part.m_softening = std::stod(word);
            } else if (i == 9) {
                part.m_potential = std::stod(word);
            }
            i++;
        }
        part.update_origin_dist();
        particles.push_back(part);

        // in case the file is formatted in an invalid way
        if (i != 10) {
            return std::nullopt;
        }
        part.update_origin_dist();
        particles.emplace_back(part);
    }

    file.close();
    return particles;
}
