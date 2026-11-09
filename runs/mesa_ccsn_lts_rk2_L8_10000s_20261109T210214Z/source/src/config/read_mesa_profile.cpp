/**
 * @file read_mesa_profile.cpp
 * @brief Reconstruct radial faces from MESA shell data and extend or truncate the mesh.
 */
#include "read_mesa_profile.h"
#include "../physcon/unit.h"
#include "../utils/utils.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace config
{
namespace
{
    /// Largest active cell count representable by MeshBlock.
    constexpr std::size_t max_cells = 10000;

    /**
     * @brief Compute the width of a radial interval in natural-log space.
     * @param[in] lo Positive inner radius.
     * @param[in] hi Outer radius in the same length units as lo.
     * @pre Both radii are finite and hi is greater than lo.
     * @return The dimensionless span @f$\ln(hi/lo)@f$.
     * @throws std::invalid_argument If the inner radius is nonpositive.
     */
    double log_span(double lo, double hi)
    {
        if (lo <= 0.0)
            throw std::invalid_argument("MESA log padding requires positive radii");
        const double width = hi - lo; // Radial width; width / lo is the relative increase in radius.
        if (width < lo)
            return std::log1p(width / lo); // Preserve small relative spacings in thin shells.
        return std::log(hi) - std::log(lo);
    }

    /**
     * @brief Determine the padding cell count on one side of the MESA grid.
     * @param[in] lo Inner face of the padding interval in FLARE length units.
     * @param[in] hi Outer face of the interval. Equal endpoints require zero cells.
     * @param[in] ref_lo Inner face of the adjacent native shell, used to estimate spacing.
     * @param[in] ref_hi Outer face of that shell, in the same units as lo.
     * @param[in] logarithmic Selects log-radius spans when true and radial widths when false.
     * @param[in] requested Explicit cell count. Zero enables estimation from the reference shell.
     * @return Zero for an empty interval, otherwise the requested count or the rounded
     *         ratio of interval width to reference width, with a minimum of one cell.
     * @throws std::invalid_argument If a logarithmic span needs a nonpositive radius.
     * @throws std::length_error If the estimated count exceeds the mesh index range.
     */
    std::size_t padding_count(double lo, double hi, double ref_lo, double ref_hi,
                              bool logarithmic, std::size_t requested)
    {
        if (lo == hi)
            return 0;
        const double span = logarithmic ? log_span(lo, hi) : hi - lo; // Padding interval width on the chosen scale.
        if (requested != 0)
            return requested;
        if (logarithmic && ref_lo == 0.0)
            throw std::invalid_argument("Specify a MESA padding count when the reference zone starts at zero");
        const double reference = logarithmic ? log_span(ref_lo, ref_hi) : ref_hi - ref_lo; // Adjacent shell width on the same scale.
        const double count = std::round(span / reference); // Check the rounded estimate before converting to an integer.
        if (!std::isfinite(count) || count > max_cells)
            throw std::length_error("Too many MESA padding cells; specify a smaller count");
        return static_cast<std::size_t>(std::max(1.0, count));
    }

    /**
     * @brief Append padding faces from the last stored face to the requested endpoint.
     * @param[in,out] faces Increasing radii in FLARE length units. The last stored radius is
     *                the starting face. The function appends one face per added cell.
     * @param[in] end Outer radius of the padding interval. Used exactly for the last face.
     * @param[in] count Number of cells to append. Zero leaves the face list unchanged.
     * @param[in] logarithmic Selects equal radius ratios when true and equal widths when false.
     * @pre For a positive count, faces is nonempty and end exceeds its last element.
     * @throws std::invalid_argument If logarithmic spacing starts at zero or below, or
     *         adjacent faces coincide at working precision.
     */
    void append_padding(std::vector<double>& faces, double end, std::size_t count,
                        bool logarithmic)
    {
        if (count == 0)
            return;
        const double start = faces.back(); // Left endpoint already in the grid.
        const double span = logarithmic ? log_span(start, end) : end - start; // Total interval on the selected scale.
        // i is the local index of an interior face, measured from the start.
        for (std::size_t i = 1; i < count; ++i)
        {
            const double t = static_cast<double>(i) / static_cast<double>(count); // Fractional position on the selected scale.
            double radius; // Face coordinate in FLARE length units.
            if (!logarithmic)
                radius = std::fma(t, span, start);
            else if (span < 1.0)
                radius = std::fma(start, std::expm1(t * span), start);
            else
                radius = std::exp(std::fma(t, span, std::log(start)));
            if (!(radius > faces.back() && radius < end))
                throw std::invalid_argument("MESA padding faces coincide at working precision");
            faces.push_back(radius);
        }
        faces.push_back(end);
    }
}

MesaProfile read_mesa_profile(const std::string& filename, double radius_inner,
                              double radius_outer, const std::string& radius_scale,
                              std::size_t inner_cells, std::size_t outer_cells,
                              double solar_radius_cm)
{
    // The scale controls only the padding intervals. Native MESA faces keep their positions.
    const bool logarithmic = utils::iequals(radius_scale, "log");
    if (!logarithmic && !utils::iequals(radius_scale, "linear"))
        throw std::invalid_argument("MESA radius_scale must be linear or log");
    if (!std::isfinite(radius_inner) || !std::isfinite(radius_outer))
        throw std::invalid_argument("MESA radius bounds must be finite");
    // Number of centimeters in one FLARE length unit.
    const double unit_length = physcon::units.unit_length;
    if (!std::isfinite(unit_length) || unit_length <= 0.0 ||
        !std::isfinite(solar_radius_cm) || solar_radius_cm <= 0.0)
        throw std::invalid_argument("MESA length units must be finite and positive");

    // Full parameter rows are reversed together so their columns stay aligned.
    Table zones = read_table(filename, 5, true);
    if (zones.data.empty())
        throw std::invalid_argument("MESA profile contains no zones: " + filename);
    std::reverse(zones.data.begin(), zones.data.end()); // MESA writes surface to center.
    // Each r starts as log10(R/Rsun) and becomes a native outer-face radius in code units.
    auto outer_faces = zones.column("logR");
    for (double& r : outer_faces)
        r = (std::pow(10.0, r) * solar_radius_cm) / unit_length;

    // logR is the OUTER face of each zone. Recover the one missing face from dr (cm).
    const double first_outer = outer_faces.front();
    const double first_width = zones.column("dr").front() / unit_length;
    if (!std::isfinite(first_width) || first_width <= 0.0)
        throw std::invalid_argument("Invalid innermost MESA zone width: " + filename);
    double mesa_inner = first_outer - first_width;
    // Allow roundoff when subtracting first_width from first_outer near the origin.
    const double tolerance = 16 * std::numeric_limits<double>::epsilon() * first_outer;
    if (mesa_inner < -tolerance)
        throw std::invalid_argument("Negative MESA inner face; check logR, dr and solar_radius_cm");
    mesa_inner = std::max(0.0, mesa_inner);
    double previous = mesa_inner; // Last radius accepted in the monotonicity check.
    for (double r : outer_faces)
    {
        if (!std::isfinite(r) || r <= previous)
            throw std::invalid_argument("MESA shell faces must increase from center to surface: " + filename);
        previous = r;
    }

    // mesa_inner and mesa_outer bound the native model. inner and outer are the resolved targets.
    const double mesa_outer = outer_faces.back();
    const double inner = radius_inner < 0.0 ? mesa_inner : radius_inner;
    const double outer = radius_outer < 0.0 ? mesa_outer : radius_outer;
    if (inner > mesa_outer || outer < mesa_inner || outer <= inner)
        throw std::invalid_argument("Require inner < outer, inner <= MESA surface and outer >= MESA inner face");

    // Keep shells with positive radial overlap. At an exact inner face the shell
    // ending there is excised; at an exact outer face the shell ending there is kept.
    const double profile_begin = std::max(inner, mesa_inner);
    const double profile_end = std::min(outer, mesa_outer);
    const auto first_shell = std::upper_bound(outer_faces.begin(), outer_faces.end(), profile_begin);
    const auto cut = std::lower_bound(first_shell, outer_faces.end(), profile_end);
    std::size_t kept = 0;
    if (profile_end > profile_begin)
        kept = static_cast<std::size_t>(cut - first_shell) + 1;

    // last_inner bounds the outermost native shell from below for outer spacing estimates.
    const double last_inner = outer_faces.size() == 1 ? mesa_inner : outer_faces[outer_faces.size() - 2];
    inner_cells = padding_count(inner, profile_begin, mesa_inner, first_outer, logarithmic, inner_cells);
    outer_cells = padding_count(profile_end, outer, last_inner, mesa_outer, logarithmic, outer_cells);
    if (kept > max_cells || inner_cells > max_cells - kept || outer_cells > max_cells - kept - inner_cells)
        throw std::length_error("MESA radial grid exceeds the MeshBlock cell count range");

    // Assemble active faces in order: inner padding, retained MESA shells, outer padding.
    std::vector<double> faces;
    faces.reserve(inner_cells + kept + outer_cells + 1);
    faces.push_back(inner);
    append_padding(faces, profile_begin, inner_cells, logarithmic);
    if (kept > 0)
    {
        faces.insert(faces.end(), first_shell, cut);
        faces.push_back(profile_end);
    }
    append_padding(faces, outer, outer_cells, logarithmic);

    // Apply the same inner cut to the parameter rows as to the native shells.
    zones.data.erase(zones.data.begin(), zones.data.begin() + (first_shell - outer_faces.begin()));
    if (inner_cells > 0)
    {
        // Save the source parameters before resizing, including when the cut retains no native shell.
        const auto innermost_row = zones.data.front();
        zones.data.resize(kept);
        zones.data.insert(zones.data.begin(), inner_cells, innermost_row);
    }
    else
        zones.data.resize(kept);
    return {std::move(faces), std::move(zones)};
}
}
