#include "data.hpp"
#include "logging.hpp"

auto main(int argc, char *argv[]) -> int {
    if (argc != 2) {
        Logging::err("File argument missing");
        return -1;
    }

    auto particles_opt = Data::read_data(argv[1]);
    if (not particles_opt.has_value()) {
        Logging::err("Error while reading file");
        return -1;
    };

    auto particles = particles_opt.value();
}
