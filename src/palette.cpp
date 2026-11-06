#include "mandelbrot/palette.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace mandelbrot {

namespace {

struct Stop {
    double position;
    Rgb color;
};

// Each gradient starts and ends with the same color so it can be repeated seamlessly.
constexpr std::array<Stop, 6> kClassic{{
    {0.0000, {0, 7, 100}},
    {0.1600, {32, 107, 203}},
    {0.4200, {237, 255, 255}},
    {0.6425, {255, 170, 0}},
    {0.8575, {0, 2, 0}},
    {1.0000, {0, 7, 100}},
}};

constexpr std::array<Stop, 6> kFire{{
    {0.00, {20, 0, 0}},
    {0.25, {160, 20, 0}},
    {0.50, {255, 140, 0}},
    {0.70, {255, 240, 160}},
    {0.85, {120, 30, 0}},
    {1.00, {20, 0, 0}},
}};

constexpr std::array<Stop, 6> kOcean{{
    {0.00, {0, 10, 30}},
    {0.30, {0, 80, 140}},
    {0.55, {40, 200, 210}},
    {0.70, {230, 255, 250}},
    {0.85, {0, 60, 110}},
    {1.00, {0, 10, 30}},
}};

constexpr std::array<Stop, 3> kGrayscale{{
    {0.0, {0, 0, 0}},
    {0.5, {255, 255, 255}},
    {1.0, {0, 0, 0}},
}};

std::uint8_t lerp(std::uint8_t a, std::uint8_t b, double t) {
    return static_cast<std::uint8_t>(std::lround(a + (b - a) * t));
}

template <std::size_t N> Rgb sample(const std::array<Stop, N>& stops, double position) {
    for (std::size_t i = 1; i < N; ++i) {
        if (position <= stops[i].position) {
            const Stop& lo = stops[i - 1];
            const Stop& hi = stops[i];
            const double t = (position - lo.position) / (hi.position - lo.position);
            return {lerp(lo.color.r, hi.color.r, t), lerp(lo.color.g, hi.color.g, t),
                    lerp(lo.color.b, hi.color.b, t)};
        }
    }
    return stops[N - 1].color;
}

} // namespace

std::optional<Palette> parse_palette(std::string_view name) {
    if (name == "classic") {
        return Palette::Classic;
    }
    if (name == "fire") {
        return Palette::Fire;
    }
    if (name == "ocean") {
        return Palette::Ocean;
    }
    if (name == "grayscale" || name == "greyscale") {
        return Palette::Grayscale;
    }
    return std::nullopt;
}

std::string_view palette_name(Palette palette) {
    switch (palette) {
    case Palette::Classic:
        return "classic";
    case Palette::Fire:
        return "fire";
    case Palette::Ocean:
        return "ocean";
    case Palette::Grayscale:
        return "grayscale";
    }
    return "classic";
}

Rgb palette_color(Palette palette, double smooth_iterations) {
    // A square-root mapping keeps the color bands a similar width at both shallow and deep zooms.
    const double cycle = std::sqrt(std::max(0.0, smooth_iterations)) / 4.0;
    const double position = cycle - std::floor(cycle);

    switch (palette) {
    case Palette::Classic:
        return sample(kClassic, position);
    case Palette::Fire:
        return sample(kFire, position);
    case Palette::Ocean:
        return sample(kOcean, position);
    case Palette::Grayscale:
        return sample(kGrayscale, position);
    }
    return sample(kClassic, position);
}

} // namespace mandelbrot
