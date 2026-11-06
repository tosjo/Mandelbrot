#pragma once

#include "mandelbrot/image.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace mandelbrot {

/// CRC-32 as used by PNG and zlib.
std::uint32_t crc32(const std::uint8_t* data, std::size_t size, std::uint32_t crc = 0);

/// Adler-32 checksum as used by the zlib container format.
std::uint32_t adler32(const std::uint8_t* data, std::size_t size);

/// Compresses data into a zlib stream (DEFLATE with LZ77 and fixed Huffman codes).
std::vector<std::uint8_t> zlib_compress(const std::vector<std::uint8_t>& data);

/// Encodes an image as an 8-bit truecolor PNG file.
std::vector<std::uint8_t> encode_png(const Image& image);

/// Encodes an image as a binary PPM (P6) file.
std::vector<std::uint8_t> encode_ppm(const Image& image);

/// Writes the image to path, choosing PNG or PPM from the file extension (PNG by default).
/// Returns false if the file could not be written.
bool write_image(const std::string& path, const Image& image);

} // namespace mandelbrot
