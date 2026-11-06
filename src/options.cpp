#include "mandelbrot/options.hpp"

#include <charconv>
#include <cmath>
#include <cstdlib>
#include <optional>
#include <string>

namespace mandelbrot {

namespace {

template <typename T> std::optional<T> parse_integer(std::string_view text, T min, T max) {
    T value{};
    const char* end = text.data() + text.size();
    const auto [ptr, ec] = std::from_chars(text.data(), end, value);
    if (ec != std::errc{} || ptr != end || value < min || value > max) {
        return std::nullopt;
    }
    return value;
}

std::optional<double> parse_double(std::string_view text) {
    if (text.empty()) {
        return std::nullopt;
    }
    const std::string copy(text);
    char* end = nullptr;
    const double value = std::strtod(copy.c_str(), &end);
    if (end != copy.c_str() + copy.size() || !std::isfinite(value)) {
        return std::nullopt;
    }
    return value;
}

// Splits "a<sep>b" into its two halves.
std::optional<std::pair<std::string_view, std::string_view>> split_pair(std::string_view text,
                                                                        char separator) {
    const auto pos = text.find(separator);
    if (pos == std::string_view::npos) {
        return std::nullopt;
    }
    return std::pair{text.substr(0, pos), text.substr(pos + 1)};
}

constexpr std::size_t kMaxDimension = 32768;

} // namespace

ParseResult parse_arguments(const std::vector<std::string_view>& args) {
    ParseResult result;
    Options& options = result.options;
    RenderSettings& settings = options.settings;

    auto fail = [&](std::string message) {
        result.error = std::move(message);
        return result;
    };

    for (std::size_t i = 0; i < args.size(); ++i) {
        std::string_view arg = args[i];
        std::optional<std::string_view> inline_value;
        if (arg.size() > 2 && arg.substr(0, 2) == "--") {
            if (const auto eq = arg.find('='); eq != std::string_view::npos) {
                inline_value = arg.substr(eq + 1);
                arg = arg.substr(0, eq);
            }
        }

        auto is = [&](std::string_view short_name, std::string_view long_name) {
            return (!short_name.empty() && arg == short_name) || arg == long_name;
        };
        auto value = [&]() -> std::optional<std::string_view> {
            if (inline_value) {
                return inline_value;
            }
            if (i + 1 < args.size()) {
                return args[++i];
            }
            return std::nullopt;
        };
        auto flag = [&]() { return !inline_value; };

        if (is("-h", "--help") && flag()) {
            options.show_help = true;
        } else if (is("", "--version") && flag()) {
            options.show_version = true;
        } else if (is("-q", "--quiet") && flag()) {
            options.quiet = true;
        } else if (is("-s", "--size")) {
            const auto text = value();
            const auto parts = text ? split_pair(*text, 'x') : std::nullopt;
            const auto width =
                parts ? parse_integer<std::size_t>(parts->first, 1, kMaxDimension) : std::nullopt;
            const auto height =
                parts ? parse_integer<std::size_t>(parts->second, 1, kMaxDimension) : std::nullopt;
            if (!width || !height) {
                return fail("--size expects WIDTHxHEIGHT, e.g. 1920x1080 (max " +
                            std::to_string(kMaxDimension) + " per side)");
            }
            settings.width = *width;
            settings.height = *height;
        } else if (is("-c", "--center")) {
            const auto text = value();
            const auto parts = text ? split_pair(*text, ',') : std::nullopt;
            const auto re = parts ? parse_double(parts->first) : std::nullopt;
            const auto im = parts ? parse_double(parts->second) : std::nullopt;
            if (!re || !im) {
                return fail("--center expects RE,IM, e.g. -0.75,0.1");
            }
            settings.center_re = *re;
            settings.center_im = *im;
        } else if (is("-z", "--zoom")) {
            const auto text = value();
            const auto zoom = text ? parse_double(*text) : std::nullopt;
            if (!zoom || *zoom <= 0.0) {
                return fail("--zoom expects a positive number");
            }
            settings.zoom = *zoom;
        } else if (is("-i", "--iterations")) {
            const auto text = value();
            const auto iterations =
                text ? parse_integer<unsigned>(*text, 1, 100'000'000) : std::nullopt;
            if (!iterations) {
                return fail("--iterations expects a number between 1 and 100000000");
            }
            settings.max_iterations = *iterations;
        } else if (is("-t", "--threads")) {
            const auto text = value();
            const auto threads = text ? parse_integer<unsigned>(*text, 0, 1024) : std::nullopt;
            if (!threads) {
                return fail("--threads expects a number between 0 (auto) and 1024");
            }
            settings.threads = *threads;
        } else if (is("-a", "--antialias")) {
            const auto text = value();
            const auto samples = text ? parse_integer<unsigned>(*text, 1, 8) : std::nullopt;
            if (!samples) {
                return fail("--antialias expects a number between 1 and 8");
            }
            settings.antialias = *samples;
        } else if (is("-p", "--palette")) {
            const auto text = value();
            const auto palette = text ? parse_palette(*text) : std::nullopt;
            if (!palette) {
                return fail("--palette expects one of: classic, fire, ocean, grayscale");
            }
            settings.palette = *palette;
        } else if (is("-o", "--output")) {
            const auto text = value();
            if (!text || text->empty()) {
                return fail("--output expects a file name");
            }
            options.output = std::string(*text);
        } else {
            return fail("unknown option '" + std::string(args[i]) + "' (see --help)");
        }
    }
    return result;
}

std::string usage(std::string_view program_name) {
    std::string text = "Usage: ";
    text += program_name;
    text += R"( [options]

Renders the Mandelbrot set to a PNG (or PPM) image.

Options:
  -s, --size WxH         image size in pixels                 (default 1200x800)
  -c, --center RE,IM     center of the view                   (default -0.75,0)
  -z, --zoom FACTOR      magnification; 1 shows the whole set (default 1)
  -i, --iterations N     maximum iterations per point         (default 1000)
  -a, --antialias N      supersample NxN per pixel, 1-8       (default 1)
  -p, --palette NAME     classic, fire, ocean or grayscale    (default classic)
  -t, --threads N        worker threads, 0 = all cores        (default 0)
  -o, --output FILE      output file, .png or .ppm            (default mandelbrot.png)
  -q, --quiet            don't print progress information
  -h, --help             show this help and exit
      --version          show the version and exit

Example:
  )";
    text += program_name;
    text += " -s 1920x1080 -c -0.743643887,0.131825904 -z 5000 -i 5000 -a 2\n";
    return text;
}

} // namespace mandelbrot
