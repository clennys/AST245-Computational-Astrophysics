#ifndef COMPASTRO_SYSTEM_H_
#define COMPASTRO_SYSTEM_H_

#include "particle.hpp"
#include "shell.hpp"

#include <limits>

class System {
  public:
    //=============================================================================================
    // Static variables and constants
    //=============================================================================================
    static constexpr double k_mean_inter_dist = 3.6078845976382845;
    /// Non dimensional mass of one
    static constexpr double k_non_dim_mass = 1.;
    // static constexpr double km_non_dim_mass = 92.4259;
    /// Softening to be applied in Force calculcation
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

    explicit System(const std::string_view &path_name);

    System();
    System(System &&) = default;
    System(const System &) = default;
    System &operator=(System &&) = default;
    System &operator=(const System &) = default;
    ~System();

    auto precalc_mean_inter_part_dist() -> double;
    auto precalc_consts() -> void;

    auto transform_vectors()
        -> std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>;

    /// @brief Return the particle that is the furthest away
    /// @note In a running system, the distances need to be calculated at each step
    auto get_max_distance() -> Particle3D;

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
    auto density_hernquist(const double rad) const -> double;
    auto newton_force(const double rad) const -> double;

    /// Return the mass found within a radius, not using `Histogram` or `Shells`
    auto get_constrained_shell_mass(const double lower_rad, const double upper_rad) const -> double;

    auto calc_direct_force() -> void;

    /// Helper method to adjust a radius to a bin size
    auto convert_lin_to_log(const int no_bins, const double val) const -> double;
    /// Helper method to add a minimal epsilon to values to circumvent log(0) errors
    static auto fit_log_to_plot(const double val) -> double;

  private:
};

#endif // ! COMPASTRO_SYSTEM_H_
