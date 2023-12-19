#include <iostream>

#include "data.hpp"

auto main(int argc, char *argv[]) -> int {
    if (!argc) {
        std::cerr << "Needed argument\n";
        return -1;
    }

    auto particles_opt = Data::read_data(argv[1]);
    if (not particles_opt.has_value()) {
        std::cerr << "Error while reading file\n";
        return -1;
    };
    auto particles = particles_opt.value();
}
