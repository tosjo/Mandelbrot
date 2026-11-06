#pragma once

#include "mandelbrot/image.hpp"

#include <optional>
#include <string_view>

namespace mandelbrot {

enum class Palette { Classic, Fire, Ocean, Grayscale };

/// Parses a palette name ("classic", "fire", "ocean", "grayscale"). Case-sensitive.
std::optional<Palette> parse_palette(std::string_view name);

/// Returns the canonical name of a palette.
std::string_view palette_name(Palette palette);

/// Maps a smooth (fractional) escape count to a color. The mapping is cyclic, so arbitrarily
/// deep zooms keep cycling through the gradient instead of saturating.
Rgb palette_color(Palette palette, double smooth_iterations);

} // namespace mandelbrot
