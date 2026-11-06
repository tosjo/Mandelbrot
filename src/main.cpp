#include "mandelbrot/options.hpp"
#include "mandelbrot/png.hpp"
#include "mandelbrot/render.hpp"

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string_view>
#include <vector>

#ifndef MANDELBROT_VERSION
#define MANDELBROT_VERSION "2.0.0"
#endif

int main(int argc, char* argv[]) {
    const std::string_view program = argc > 0 ? argv[0] : "mandelbrot";
    std::vector<std::string_view> args(argv + (argc > 0 ? 1 : 0), argv + argc);

    const mandelbrot::ParseResult parsed = mandelbrot::parse_arguments(args);
    if (!parsed.ok()) {
        std::cerr << "error: " << parsed.error << '\n';
        return EXIT_FAILURE;
    }
    const mandelbrot::Options& options = parsed.options;
    if (options.show_help) {
        std::cout << mandelbrot::usage(program);
        return EXIT_SUCCESS;
    }
    if (options.show_version) {
        std::cout << "mandelbrot " << MANDELBROT_VERSION << '\n';
        return EXIT_SUCCESS;
    }

    const auto start = std::chrono::steady_clock::now();
    const mandelbrot::Image image = mandelbrot::render(options.settings);
    const auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start);

    if (!mandelbrot::write_image(options.output, image)) {
        std::cerr << "error: could not write '" << options.output << "'\n";
        return EXIT_FAILURE;
    }

    if (!options.quiet) {
        std::cout << "Rendered " << image.width() << 'x' << image.height() << " in "
                  << elapsed.count() << " s -> " << options.output << '\n';
    }
    return EXIT_SUCCESS;
}
