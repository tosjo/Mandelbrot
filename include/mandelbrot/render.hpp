#pragma once

#include "mandelbrot/image.hpp"
#include "mandelbrot/palette.hpp"

#include <cstddef>

namespace mandelbrot {

/// Value returned by escape_time() for points that never escaped (i.e. are in the set).
inline constexpr double kInside = -1.0;

/// Computes the smooth (normalized) iteration count for the point c = re + im*i.
/// Returns kInside if the orbit did not escape within max_iterations, otherwise a value >= 0.
double escape_time(double re, double im, unsigned max_iterations);

struct RenderSettings {
    std::size_t width = 1200;
    std::size_t height = 800;
    double center_re = -0.75;
    double center_im = 0.0;
    /// Zoom factor; at zoom 1 the shorter side of the image spans 2.5 units of the complex plane.
    double zoom = 1.0;
    unsigned max_iterations = 1000;
    /// Number of worker threads; 0 means "use all hardware threads".
    unsigned threads = 0;
    /// Supersampling factor per axis (1 = no anti-aliasing, 3 = 9 samples per pixel).
    unsigned antialias = 1;
    Palette palette = Palette::Classic;
};

/// Distance between neighbouring pixel centers in the complex plane.
double pixel_size(const RenderSettings& settings);

/// Renders the Mandelbrot set using the given settings. Rows are distributed dynamically
/// over the worker threads, so expensive rows near the set don't stall the others.
Image render(const RenderSettings& settings);

} // namespace mandelbrot
