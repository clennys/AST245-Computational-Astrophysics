#ifndef COMPASTRO_SYSTEM_H_
#define COMPASTRO_SYSTEM_H_

#include "node.hpp"
#include "particle.hpp"
#include "shell.hpp"

#include <limits>

class System {
  public:
    //=============================================================================================
    // Static variables and constants
    //=============================================================================================

    /// Graviational Constant in pc * (km/s)^2 / M_\odot
    static constexpr double k_G = 4.3009172706e-3;
    /// Precalculated Mean inter-particle distance
    static constexpr double k_mean_inter_dist = 4.023775510528517;
    /// Non dimensional particle mass of one
    static constexpr double k_non_dim_mass = 1.;
    /// Dimensional Mass per particle
    static constexpr double k_dim_mass = 92.4259;
    /// Array of divisors to use for softening
    static constexpr auto softening_divisors = {1., 10., 20., 50., 100., 200.};
    /// Softening to be applied in Force calculcation. Made static, in order for other Classes to
    /// access it
    static double s_softening_length;
    static double s_softening;

    //=============================================================================================
    // Regular member variables
    //=============================================================================================
    PartVec m_particles;
    double m_total_mass = 0.;
    double m_half_mass_rad = 0.;
    double m_scale_length = 0.;
    double m_min_rad = std::numeric_limits<double>::max();
    double m_max_rad = 0.;
    double m_softening = 0.;
    double m_relaxation = 0.;

    //=============================================================================================
    // Ctors, Dtors, etc.
    //=============================================================================================
    System() = default;
    System(System &&) = default;
    System(const System &) = default;
    System &operator=(System &&) = default;
    System &operator=(const System &) = default;
    ~System() = default;

    /// Special Init method, initializes particles with `path_name` file and calculates constant for
    auto init_system(const std::string_view &path_name) -> void;

    /// Return the count of particles in the whole system cast as an `int` for correct calculations
    /// with it
    auto system_int_size() const -> int;

    /// Calculate the Mean inter-particle distance
    ///
    /// NOTE: Recommended to be run once, as O(N^2)
    [[nodiscard]] auto precalc_mean_inter_part_dist() -> double;

    auto transform_vectors()
        -> std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>;

    /// @brief Return the particle that is the furthest away
    /// @note In a running system, the distances need to be calculated at each step
    [[nodiscard]] auto get_max_distance() -> Particle3D;

    /// Calculate and return the total mass inside the system
    [[nodiscard]] auto calc_total_mass() const -> double;
    /// Calculate and set the total mass variable
    auto update_total_mass() -> void;

    /// Calculate the half mass and return it
    [[nodiscard]] auto calc_half_mass_radius(const ShellVec &shells) const -> double;
    /// Calculate the half mass and set the half mass member value
    auto update_half_mass_radius(const ShellVec &shells) -> void;

    /// Calculate the scale length and return it
    [[nodiscard]] auto calc_scale_length() const -> double;
    /// Calculate the scale length and set its member variable
    auto update_scale_length() -> void;

    /// update the minimal radius by comparing it the the one passed to it
    auto update_min_rad(const double rad) -> void;
    /// update the maximal radius by comparing it the the one passed to it
    auto update_max_rad(const double rad) -> void;

    /// return the analytical density profile within a radius for Hernquist
    [[nodiscard]] auto density_hernquist(const double rad) const -> double;
    [[nodiscard]] auto newton_force(const double rad) const -> double;

    auto calc_direct_initial_force() -> void;

    /// Helper method to adjust a radius to a bin size
    ///
    /// The ratio is created by `m_min_rad` and `m_max_rad`.
    /// The exponent scales the value to between 0 and 1., by dividing `val` and `no_bins`.
    /// The multiplication of the value with `m_min_rad` guarantees that the value falls in to the
    /// designated range.
    [[nodiscard]] auto convert_lin_to_log(const int no_bins, const double val) const -> double;

    /// Helper method to add a minimal epsilon to values to circumvent log(0) errors
    [[nodiscard]] static auto fit_log_to_plot(const double val) -> double;

    /// Do one step forward in the system
    auto solver_do_step(const double delta_time) -> void;

    [[nodiscard]] auto calc_relaxation() const -> double;

    auto update_relaxation() -> void;

    auto calc_overall_bounding_cube() -> BoundingCube;

  private:
    /// Calculate runtime constants for system deployment
    auto precalc_consts() -> void;
};

#endif // ! COMPASTRO_SYSTEM_H_
