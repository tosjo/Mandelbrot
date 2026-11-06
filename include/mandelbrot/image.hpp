#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace mandelbrot {

struct Rgb {
    std::uint8_t r = 0;
    std::uint8_t g = 0;
    std::uint8_t b = 0;

    friend bool operator==(const Rgb& lhs, const Rgb& rhs) {
        return lhs.r == rhs.r && lhs.g == rhs.g && lhs.b == rhs.b;
    }
    friend bool operator!=(const Rgb& lhs, const Rgb& rhs) { return !(lhs == rhs); }
};

/// A simple row-major 8-bit RGB image.
class Image {
  public:
    Image() = default;
    Image(std::size_t width, std::size_t height)
        : width_(width), height_(height), pixels_(width * height) {}

    std::size_t width() const { return width_; }
    std::size_t height() const { return height_; }

    Rgb& at(std::size_t x, std::size_t y) { return pixels_[y * width_ + x]; }
    const Rgb& at(std::size_t x, std::size_t y) const { return pixels_[y * width_ + x]; }

    Rgb* row(std::size_t y) { return pixels_.data() + y * width_; }
    const Rgb* row(std::size_t y) const { return pixels_.data() + y * width_; }

  private:
    std::size_t width_ = 0;
    std::size_t height_ = 0;
    std::vector<Rgb> pixels_;
};

} // namespace mandelbrot
