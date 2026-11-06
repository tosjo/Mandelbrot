#pragma once

#include "mandelbrot/render.hpp"

#include <string>
#include <string_view>
#include <vector>

namespace mandelbrot {

struct Options {
    RenderSettings settings;
    std::string output = "mandelbrot.png";
    bool quiet = false;
    bool show_help = false;
    bool show_version = false;
};

struct ParseResult {
    Options options;
    /// Empty on success, otherwise a human-readable error message.
    std::string error;

    bool ok() const { return error.empty(); }
};

/// Parses command line arguments (excluding the program name).
ParseResult parse_arguments(const std::vector<std::string_view>& args);

/// Returns the --help text.
std::string usage(std::string_view program_name);

} // namespace mandelbrot
