#include <fmt/core.h>

#include "mpint.h"

auto main() -> int {
    mpint test = 10 * (1000 * 1000);
    fmt::println("{}", test);
}
