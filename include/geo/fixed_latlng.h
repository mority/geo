#pragma once

#include <cinttypes>
#include <cmath>
#include <limits>

#include "geo/latlng.h"

namespace geo {

struct fixed_latlng {
  static constexpr auto const kCoordinatePrecision = 1'000'0000;
  static constexpr auto const kUndefined =
      std::numeric_limits<std::int32_t>::max();

  static std::int32_t double_to_fix(double const c) noexcept {
    return static_cast<std::int32_t>(std::round(c * kCoordinatePrecision));
  }
  static constexpr double fix_to_double(std::int32_t const c) noexcept {
    return static_cast<double>(c) / kCoordinatePrecision;
  }

  static fixed_latlng from_latlng(latlng const& x) {
    return {double_to_fix(x.lat()), double_to_fix(x.lng())};
  }

  double lat() const { return fix_to_double(lat_); }
  double lng() const { return fix_to_double(lng_); }

  constexpr std::int32_t x() const noexcept { return lng_; }
  constexpr std::int32_t y() const noexcept { return lat_; }

  fixed_latlng& set_x(std::int32_t const x) noexcept {
    lng_ = x;
    return *this;
  }
  fixed_latlng& set_y(std::int32_t const y) noexcept {
    lat_ = y;
    return *this;
  }

  constexpr bool valid() const noexcept {
    return lat_ >= -90 * kCoordinatePrecision &&
           lat_ <= 90 * kCoordinatePrecision &&
           lng_ >= -180 * kCoordinatePrecision &&
           lng_ <= 180 * kCoordinatePrecision;
  }

  constexpr bool is_set() const noexcept {
    return lat_ != kUndefined && lng_ != kUndefined;
  }

  constexpr bool operator==(fixed_latlng const& o) const noexcept {
    return lat_ == o.lat_ && lng_ == o.lng_;
  }
  constexpr bool operator!=(fixed_latlng const& o) const noexcept {
    return !(*this == o);
  }
  constexpr bool operator<(fixed_latlng const& o) const noexcept {
    return (lat_ == o.lat_ && lng_ < o.lng_) || lat_ < o.lat_;
  }
  constexpr bool operator>(fixed_latlng const& o) const noexcept {
    return (o.lat_ == lat_ && o.lng_ < lng_) || o.lat_ < lat_;
  }

  operator latlng() const { return {lat(), lng()}; }

  std::int32_t lat_{kUndefined};
  std::int32_t lng_{kUndefined};
};

}  // namespace geo
