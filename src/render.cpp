#include "mandelbrot/render.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <thread>
#include <vector>

namespace mandelbrot {

namespace {

// A large bailout radius makes the smooth coloring formula accurate.
constexpr double kBailoutSquared = 256.0 * 256.0;

// Points in the main cardioid or the period-2 bulb never escape; checking for them up front
// skips the most expensive part of a typical render.
bool in_main_cardioid_or_bulb(double re, double im) {
    const double im2 = im * im;
    const double x = re - 0.25;
    const double q = x * x + im2;
    if (q * (q + x) <= 0.25 * im2) {
        return true;
    }
    const double y = re + 1.0;
    return y * y + im2 <= 1.0 / 16.0;
}

unsigned resolve_thread_count(unsigned requested) {
    if (requested != 0) {
        return requested;
    }
    return std::max(1u, std::thread::hardware_concurrency());
}

} // namespace

double escape_time(double re, double im, unsigned max_iterations) {
    if (in_main_cardioid_or_bulb(re, im)) {
        return kInside;
    }

    double zr = 0.0;
    double zi = 0.0;
    double zr2 = 0.0;
    double zi2 = 0.0;
    for (unsigned n = 0; n < max_iterations; ++n) {
        zi = 2.0 * zr * zi + im;
        zr = zr2 - zi2 + re;
        zr2 = zr * zr;
        zi2 = zi * zi;
        const double magnitude2 = zr2 + zi2;
        if (magnitude2 > kBailoutSquared) {
            // Normalized iteration count: n + 1 - log2(log|z|).
            const double log_modulus = 0.5 * std::log(magnitude2);
            const double smooth = static_cast<double>(n) + 1.0 - std::log2(log_modulus);
            return std::max(0.0, smooth);
        }
    }
    return kInside;
}

double pixel_size(const RenderSettings& settings) {
    const auto shorter_side = static_cast<double>(std::min(settings.width, settings.height));
    return 2.5 / (settings.zoom * shorter_side);
}

Image render(const RenderSettings& settings) {
    Image image(settings.width, settings.height);
    if (settings.width == 0 || settings.height == 0) {
        return image;
    }

    const double step = pixel_size(settings);
    const double left = settings.center_re - 0.5 * static_cast<double>(settings.width - 1) * step;
    const double top = settings.center_im + 0.5 * static_cast<double>(settings.height - 1) * step;

    const unsigned samples = std::max(1u, settings.antialias);
    const double sub_step = step / samples;
    const double sub_offset = -0.5 * step + 0.5 * sub_step;
    const unsigned samples_per_pixel = samples * samples;

    std::atomic<std::size_t> next_row{0};

    auto worker = [&] {
        for (std::size_t y = next_row++; y < settings.height; y = next_row++) {
            Rgb* row = image.row(y);
            const double im0 = top - static_cast<double>(y) * step;
            for (std::size_t x = 0; x < settings.width; ++x) {
                const double re0 = left + static_cast<double>(x) * step;
                unsigned r = 0;
                unsigned g = 0;
                unsigned b = 0;
                for (unsigned sy = 0; sy < samples; ++sy) {
                    const double im = im0 - (sub_offset + sy * sub_step);
                    for (unsigned sx = 0; sx < samples; ++sx) {
                        const double re = re0 + sub_offset + sx * sub_step;
                        const double t = escape_time(re, im, settings.max_iterations);
                        const Rgb color = t == kInside ? Rgb{} : palette_color(settings.palette, t);
                        r += color.r;
                        g += color.g;
                        b += color.b;
                    }
                }
                const unsigned half = samples_per_pixel / 2;
                row[x] = Rgb{static_cast<std::uint8_t>((r + half) / samples_per_pixel),
                             static_cast<std::uint8_t>((g + half) / samples_per_pixel),
                             static_cast<std::uint8_t>((b + half) / samples_per_pixel)};
            }
        }
    };

    const unsigned thread_count =
        std::min<unsigned>(resolve_thread_count(settings.threads),
                           static_cast<unsigned>(std::min<std::size_t>(settings.height, 4096)));
    std::vector<std::thread> pool;
    pool.reserve(thread_count - 1);
    for (unsigned i = 1; i < thread_count; ++i) {
        pool.emplace_back(worker);
    }
    worker();
    for (auto& thread : pool) {
        thread.join();
    }
    return image;
}

} // namespace mandelbrot
