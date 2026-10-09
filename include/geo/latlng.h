#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <array>
#include <iosfwd>
#include <limits>
#include <tuple>

#include "geo/constants.h"
#include "geo/rad_deg.h"

namespace geo {

struct latlng {
  double lat() const noexcept { return lat_; }
  double lng() const noexcept { return lng_; }

  friend std::ostream& operator<<(std::ostream&, latlng const&);

  friend bool operator<(latlng const& lhs, latlng const& rhs) noexcept {
    return std::tie(lhs.lat_, lhs.lng_) < std::tie(rhs.lat_, rhs.lng_);
  }

  friend bool operator==(latlng const& lhs, latlng const& rhs) noexcept {
    auto const lat_diff = std::abs(lhs.lat_ - rhs.lat_);
    auto const lng_diff = std::abs(lhs.lng_ - rhs.lng_);
    return lat_diff < 100 * std::numeric_limits<double>::epsilon() &&
           lng_diff < 100 * std::numeric_limits<double>::epsilon();
  }

  std::array<double, 2> lnglat() const noexcept { return {lng_, lat_}; }
  std::array<float, 2> lnglat_float() const noexcept {
    return {static_cast<float>(lng_), static_cast<float>(lat_)};
  }

  double lat_{0.0}, lng_{0.0};
};

struct meters_per_degree {
  double lat_{0.0}, lng_{0.0};
};

namespace detail {

// meters per degree at latitudes -90°, -89°, ..., 90°
inline std::array<meters_per_degree, 181> make_local_radii_table() {
  constexpr auto kA = 6378137.0;  // WGS84 semi-major axis
  constexpr auto kF = 1.0 / 298.257223563;  // WGS84 flattening
  constexpr auto kE2 = kF * (2.0 - kF);  // squared eccentricity

  auto t = std::array<meters_per_degree, 181>{};
  for (auto i = std::size_t{0U}; i != t.size(); ++i) {
    auto const lat = to_rad(static_cast<double>(i) - 90.0);
    auto const s = std::sin(lat);
    auto const w = 1.0 - kE2 * s * s;
    auto const meridian = kA * (1.0 - kE2) / (w * std::sqrt(w));
    auto const prime_vertical = kA / std::sqrt(w);
    t[i] = {meridian * kPI / 180.0,
            prime_vertical * std::cos(lat) * kPI / 180.0};
  }
  return t;
}

// namespace-scope instead of function-local static: no guard check on access
inline auto const kLocalRadiiTable = make_local_radii_table();

}  // namespace detail

double distance(latlng const&, latlng const&);

inline double approx_squared_distance(
    latlng const& a, latlng const& b,
    double const approx_distance_lng_degrees) {
  auto const y = std::abs(a.lat() - b.lat()) * kApproxDistanceLatDegrees;
  auto const xdiff = std::abs(a.lng() - b.lng());
  auto const x =
      (xdiff > 180.0 ? (360.0 - xdiff) : xdiff) * approx_distance_lng_degrees;
  return x * x + y * y;
}

inline double approx_squared_distance(latlng const& a, latlng const& b,
                                      meters_per_degree const& scale) {
  auto const y = (a.lat() - b.lat()) * scale.lat_;
  auto const xdiff = std::abs(a.lng() - b.lng());
  auto const x = (xdiff > 180.0 ? (360.0 - xdiff) : xdiff) * scale.lng_;
  return x * x + y * y;
}

// WGS84 meridian (M) and prime-vertical (N cos lat) radius of curvature at the
// given latitude in meters per degree, from the local radii lookup table.
// Linear interpolation between the 1° entries:
// relative error ~h²/8 with h = 1° in radians (~0.004%)
inline meters_per_degree approx_meters_per_degree(double const lat) {
  auto const& t = detail::kLocalRadiiTable;
  auto const f = std::clamp(lat + 90.0, 0.0, 180.0);
  // std::uint32_t: double -> std::size_t has no single instruction on x86-64
  auto const i = std::min(static_cast<std::uint32_t>(f), std::uint32_t{179U});
  auto const r = f - static_cast<double>(i);
  return {t[i].lat_ + r * (t[i + 1].lat_ - t[i].lat_),
          t[i].lng_ + r * (t[i + 1].lng_ - t[i].lng_)};
}

// Equirectangular (flat-earth) distance approximation using the WGS84 meridian
// and prime-vertical radii of curvature at the mean latitude, with scale
// factors from a linearly interpolated 1° lookup table.
// Error grows quadratically with the distance between the points, only use
// for nearby points away from the poles.
inline double approx_squared_distance_local_radii(latlng const& a,
                                                  latlng const& b) {
  return approx_squared_distance(
      a, b, approx_meters_per_degree(0.5 * (a.lat() + b.lat())));
}

double bearing(latlng const&, latlng const&);

latlng midpoint(latlng const&, latlng const&);

latlng closest_on_segment(latlng const& x, latlng const& segment_from,
                          latlng const& segment_to);

std::pair<latlng, double> approx_closest_on_segment(
    latlng const& x, latlng const& segment_from, latlng const& segment_to,
    double approx_distance_lng_degrees);

double lower_bound_distance_lng_degrees(latlng const&);

double approx_distance_lng_degrees(latlng const&);

latlng destination_point(latlng const& source, double const distance,
                         double const bearing);

uint32_t tile_hash_32(latlng const&);

std::uint64_t morton_encode(latlng const&);

}  // namespace geo

#if __has_include("fmt/format.h")

#include "fmt/format.h"

template <>
struct fmt::formatter<geo::latlng> : nested_formatter<double> {
  auto format(geo::latlng const& p, format_context& ctx) const {
    return write_padded(ctx, [&](auto out) {
      return format_to(out, "({}, {})", nested(p.lat_), nested(p.lng_));
    });
  }
};

#endif