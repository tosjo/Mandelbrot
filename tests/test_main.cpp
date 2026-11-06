// Small self-contained test runner, so the project has no external dependencies.

#include "mandelbrot/options.hpp"
#include "mandelbrot/palette.hpp"
#include "mandelbrot/png.hpp"
#include "mandelbrot/render.hpp"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

namespace {

int g_failures = 0;

#define CHECK(condition)                                                                           \
    do {                                                                                           \
        if (!(condition)) {                                                                        \
            std::cerr << __FILE__ << ':' << __LINE__ << ": check failed: " #condition "\n";        \
            ++g_failures;                                                                          \
        }                                                                                          \
    } while (false)

using namespace mandelbrot;

std::vector<std::uint8_t> bytes(std::string_view text) {
    return {text.begin(), text.end()};
}

void test_points_inside_the_set() {
    CHECK(escape_time(0.0, 0.0, 1000) == kInside);
    CHECK(escape_time(-1.0, 0.0, 1000) == kInside);
    CHECK(escape_time(-0.1, 0.1, 1000) == kInside);
    CHECK(escape_time(-1.75, 0.0, 1000) == kInside); // period-3 bulb on the real axis
    CHECK(escape_time(0.25, 0.0, 1000) == kInside);  // cusp of the main cardioid
    CHECK(escape_time(-2.0, 0.0, 1000) == kInside);  // tip of the antenna
}

void test_points_outside_the_set() {
    CHECK(escape_time(1.0, 0.0, 1000) >= 0.0);
    CHECK(escape_time(0.0, 2.0, 1000) >= 0.0);
    CHECK(escape_time(-2.1, 0.0, 1000) >= 0.0);
    CHECK(escape_time(0.26, 0.0, 10000) >= 0.0);
    // Points further away escape sooner.
    CHECK(escape_time(10.0, 10.0, 1000) < escape_time(0.5, 0.5, 1000));
}

void test_smooth_escape_is_continuous() {
    const double a = escape_time(0.4, 0.6, 1000);
    const double b = escape_time(0.4 + 1e-9, 0.6, 1000);
    CHECK(a >= 0.0 && b >= 0.0);
    CHECK(std::abs(a - b) < 1e-3);
}

void test_render_dimensions_and_symmetry() {
    RenderSettings settings;
    settings.width = 41;
    settings.height = 31;
    settings.center_im = 0.0;
    settings.max_iterations = 200;
    settings.threads = 3;
    const Image image = render(settings);
    CHECK(image.width() == 41);
    CHECK(image.height() == 31);
    // The set is symmetric around the real axis and the view is centered on it.
    bool symmetric = true;
    for (std::size_t y = 0; y < image.height(); ++y) {
        for (std::size_t x = 0; x < image.width(); ++x) {
            symmetric = symmetric && image.at(x, y) == image.at(x, image.height() - 1 - y);
        }
    }
    CHECK(symmetric);
    // The center pixel (-0.75, 0) is inside the set and rendered black.
    CHECK(image.at(20, 15) == Rgb{});
}

void test_render_is_independent_of_thread_count() {
    RenderSettings settings;
    settings.width = 64;
    settings.height = 40;
    settings.max_iterations = 300;
    settings.antialias = 2;
    settings.threads = 1;
    const Image single = render(settings);
    settings.threads = 7;
    const Image multi = render(settings);
    bool equal = true;
    for (std::size_t y = 0; y < single.height(); ++y) {
        for (std::size_t x = 0; x < single.width(); ++x) {
            equal = equal && single.at(x, y) == multi.at(x, y);
        }
    }
    CHECK(equal);
}

void test_palettes() {
    CHECK(parse_palette("classic") == Palette::Classic);
    CHECK(parse_palette("fire") == Palette::Fire);
    CHECK(parse_palette("ocean") == Palette::Ocean);
    CHECK(parse_palette("grayscale") == Palette::Grayscale);
    CHECK(!parse_palette("rainbow"));
    for (Palette p : {Palette::Classic, Palette::Fire, Palette::Ocean, Palette::Grayscale}) {
        CHECK(parse_palette(palette_name(p)) == p);
    }
    // Grayscale is black at the start of its cycle and white halfway.
    CHECK(palette_color(Palette::Grayscale, 0.0) == Rgb{});
    CHECK(palette_color(Palette::Grayscale, 4.0) == (Rgb{255, 255, 255}));
}

void test_checksums() {
    const auto check = bytes("123456789");
    CHECK(crc32(check.data(), check.size()) == 0xCBF43926u);
    const auto wiki = bytes("Wikipedia");
    CHECK(adler32(wiki.data(), wiki.size()) == 0x11E60398u);
}

void test_zlib_stream_layout() {
    std::vector<std::uint8_t> data(10000, 'a');
    const auto compressed = zlib_compress(data);
    CHECK(compressed.size() > 6);
    CHECK(compressed.size() < 200); // long runs must be compressed by back-references
    CHECK((compressed[0] * 256 + compressed[1]) % 31 == 0);
    const std::size_t n = compressed.size();
    const std::uint32_t trailer = (std::uint32_t{compressed[n - 4]} << 24) |
                                  (std::uint32_t{compressed[n - 3]} << 16) |
                                  (std::uint32_t{compressed[n - 2]} << 8) | compressed[n - 1];
    CHECK(trailer == adler32(data.data(), data.size()));
}

void test_png_structure() {
    Image image(3, 2);
    image.at(1, 1) = Rgb{255, 0, 0};
    const auto png = encode_png(image);
    const std::vector<std::uint8_t> signature{0x89, 'P', 'N', 'G', '\r', '\n', 0x1A, '\n'};
    CHECK(std::equal(signature.begin(), signature.end(), png.begin()));
    CHECK(std::string(png.begin() + 12, png.begin() + 16) == "IHDR");
    CHECK(png[19] == 3 && png[23] == 2); // width and height (big-endian)
    CHECK(std::string(png.end() - 8, png.end() - 4) == "IEND");
}

void test_ppm_structure() {
    Image image(2, 1);
    image.at(1, 0) = Rgb{1, 2, 3};
    const auto ppm = encode_ppm(image);
    const std::string expected_header = "P6\n2 1\n255\n";
    CHECK(std::string(ppm.begin(), ppm.begin() + 11) == expected_header);
    CHECK(ppm.size() == expected_header.size() + 6);
    CHECK(ppm[ppm.size() - 1] == 3);
}

ParseResult parse(std::vector<std::string_view> args) {
    return parse_arguments(args);
}

void test_argument_defaults() {
    const auto result = parse({});
    CHECK(result.ok());
    CHECK(result.options.output == "mandelbrot.png");
    CHECK(result.options.settings.width == 1200 && result.options.settings.height == 800);
}

void test_argument_parsing() {
    const auto result = parse({"-s", "320x200", "--center", "-0.5,0.25", "--zoom=4", "-i", "500",
                               "-t", "2", "-a", "3", "-p", "fire", "-o", "out.ppm", "-q"});
    CHECK(result.ok());
    const auto& s = result.options.settings;
    CHECK(s.width == 320 && s.height == 200);
    CHECK(s.center_re == -0.5 && s.center_im == 0.25);
    CHECK(s.zoom == 4.0);
    CHECK(s.max_iterations == 500);
    CHECK(s.threads == 2);
    CHECK(s.antialias == 3);
    CHECK(s.palette == Palette::Fire);
    CHECK(result.options.output == "out.ppm");
    CHECK(result.options.quiet);
}

void test_argument_errors() {
    CHECK(!parse({"--size", "100"}).ok());
    CHECK(!parse({"--size", "0x10"}).ok());
    CHECK(!parse({"--size", "10x-10"}).ok());
    CHECK(!parse({"--center", "1"}).ok());
    CHECK(!parse({"--center", "a,b"}).ok());
    CHECK(!parse({"--zoom", "0"}).ok());
    CHECK(!parse({"--zoom", "inf"}).ok());
    CHECK(!parse({"--iterations", "12abc"}).ok());
    CHECK(!parse({"--antialias", "9"}).ok());
    CHECK(!parse({"--palette", "rainbow"}).ok());
    CHECK(!parse({"--output"}).ok());
    CHECK(!parse({"--quiet=yes"}).ok());
    CHECK(!parse({"--frobnicate"}).ok());
}

} // namespace

int main() {
    const std::vector<std::pair<const char*, std::function<void()>>> tests{
        {"points inside the set", test_points_inside_the_set},
        {"points outside the set", test_points_outside_the_set},
        {"smooth escape is continuous", test_smooth_escape_is_continuous},
        {"render dimensions and symmetry", test_render_dimensions_and_symmetry},
        {"render is independent of thread count", test_render_is_independent_of_thread_count},
        {"palettes", test_palettes},
        {"checksums", test_checksums},
        {"zlib stream layout", test_zlib_stream_layout},
        {"png structure", test_png_structure},
        {"ppm structure", test_ppm_structure},
        {"argument defaults", test_argument_defaults},
        {"argument parsing", test_argument_parsing},
        {"argument errors", test_argument_errors},
    };

    for (const auto& [name, test] : tests) {
        const int before = g_failures;
        test();
        std::cout << (g_failures == before ? "[ OK ] " : "[FAIL] ") << name << '\n';
    }
    if (g_failures != 0) {
        std::cerr << g_failures << " check(s) failed\n";
        return EXIT_FAILURE;
    }
    std::cout << "All tests passed\n";
    return EXIT_SUCCESS;
}
