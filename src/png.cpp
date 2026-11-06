#include "mandelbrot/png.hpp"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <fstream>
#include <string_view>

namespace mandelbrot {

namespace {

using Bytes = std::vector<std::uint8_t>;

const std::array<std::uint32_t, 256>& crc_table() {
    static const auto table = [] {
        std::array<std::uint32_t, 256> t{};
        for (std::uint32_t n = 0; n < 256; ++n) {
            std::uint32_t c = n;
            for (int k = 0; k < 8; ++k) {
                c = (c & 1u) ? 0xEDB88320u ^ (c >> 1) : c >> 1;
            }
            t[n] = c;
        }
        return t;
    }();
    return table;
}

void put_u32_be(Bytes& out, std::uint32_t value) {
    out.push_back(static_cast<std::uint8_t>(value >> 24));
    out.push_back(static_cast<std::uint8_t>(value >> 16));
    out.push_back(static_cast<std::uint8_t>(value >> 8));
    out.push_back(static_cast<std::uint8_t>(value));
}

// Writes bits least-significant first, as DEFLATE requires.
class BitWriter {
  public:
    explicit BitWriter(Bytes& out) : out_(out) {}

    void write(std::uint32_t bits, unsigned count) {
        buffer_ |= static_cast<std::uint64_t>(bits) << filled_;
        filled_ += count;
        while (filled_ >= 8) {
            out_.push_back(static_cast<std::uint8_t>(buffer_));
            buffer_ >>= 8;
            filled_ -= 8;
        }
    }

    // Huffman codes are defined most-significant bit first, so they are reversed before writing.
    void write_code(std::uint32_t code, unsigned length) {
        std::uint32_t reversed = 0;
        for (unsigned i = 0; i < length; ++i) {
            reversed = (reversed << 1) | ((code >> i) & 1u);
        }
        write(reversed, length);
    }

    void flush() {
        if (filled_ > 0) {
            out_.push_back(static_cast<std::uint8_t>(buffer_));
            buffer_ = 0;
            filled_ = 0;
        }
    }

  private:
    Bytes& out_;
    std::uint64_t buffer_ = 0;
    unsigned filled_ = 0;
};

// Fixed Huffman code for a literal/length symbol (RFC 1951, section 3.2.6).
void write_literal_length_symbol(BitWriter& writer, unsigned symbol) {
    if (symbol < 144) {
        writer.write_code(0x30 + symbol, 8);
    } else if (symbol < 256) {
        writer.write_code(0x190 + (symbol - 144), 9);
    } else if (symbol < 280) {
        writer.write_code(symbol - 256, 7);
    } else {
        writer.write_code(0xC0 + (symbol - 280), 8);
    }
}

constexpr std::array<std::uint16_t, 29> kLengthBase{3,  4,  5,  6,   7,   8,   9,   10,  11, 13,
                                                    15, 17, 19, 23,  27,  31,  35,  43,  51, 59,
                                                    67, 83, 99, 115, 131, 163, 195, 227, 258};
constexpr std::array<std::uint8_t, 29> kLengthExtra{0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2,
                                                    2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0};
constexpr std::array<std::uint16_t, 30> kDistanceBase{
    1,   2,   3,   4,   5,   7,    9,    13,   17,   25,   33,   49,   65,    97,    129,
    193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577};
constexpr std::array<std::uint8_t, 30> kDistanceExtra{0, 0, 0,  0,  1,  1,  2,  2,  3,  3,
                                                      4, 4, 5,  5,  6,  6,  7,  7,  8,  8,
                                                      9, 9, 10, 10, 11, 11, 12, 12, 13, 13};

void write_match(BitWriter& writer, unsigned length, unsigned distance) {
    unsigned code = 0;
    while (code + 1 < kLengthBase.size() && kLengthBase[code + 1] <= length) {
        ++code;
    }
    write_literal_length_symbol(writer, 257 + code);
    writer.write(length - kLengthBase[code], kLengthExtra[code]);

    unsigned dcode = 0;
    while (dcode + 1 < kDistanceBase.size() && kDistanceBase[dcode + 1] <= distance) {
        ++dcode;
    }
    writer.write_code(dcode, 5);
    writer.write(distance - kDistanceBase[dcode], kDistanceExtra[dcode]);
}

void deflate_fixed(const Bytes& data, Bytes& out) {
    constexpr std::size_t kWindow = 32768;
    constexpr std::size_t kMinMatch = 3;
    constexpr std::size_t kMaxMatch = 258;
    constexpr std::size_t kMaxChain = 64;
    constexpr std::size_t kHashBits = 15;
    constexpr std::size_t kNone = static_cast<std::size_t>(-1);

    BitWriter writer(out);
    writer.write(1, 1); // BFINAL: this is the only block
    writer.write(1, 2); // BTYPE = 01: fixed Huffman codes

    std::vector<std::size_t> head(std::size_t{1} << kHashBits, kNone);
    std::vector<std::size_t> prev(kWindow, kNone);
    auto hash = [&](std::size_t i) {
        const std::uint32_t v = (std::uint32_t{data[i]} << 16) | (std::uint32_t{data[i + 1]} << 8) |
                                std::uint32_t{data[i + 2]};
        return static_cast<std::size_t>((v * 2654435761u) >> (32 - kHashBits));
    };
    auto insert = [&](std::size_t i) {
        if (i + kMinMatch <= data.size()) {
            const std::size_t h = hash(i);
            prev[i % kWindow] = head[h];
            head[h] = i;
        }
    };

    const std::size_t size = data.size();
    std::size_t i = 0;
    while (i < size) {
        std::size_t best_length = 0;
        std::size_t best_distance = 0;
        if (i + kMinMatch <= size) {
            const std::size_t max_length = std::min(kMaxMatch, size - i);
            std::size_t candidate = head[hash(i)];
            for (std::size_t chain = 0; candidate != kNone && chain < kMaxChain; ++chain) {
                const std::size_t distance = i - candidate;
                if (distance == 0 || distance > kWindow) {
                    break;
                }
                std::size_t length = 0;
                while (length < max_length && data[candidate + length] == data[i + length]) {
                    ++length;
                }
                if (length > best_length) {
                    best_length = length;
                    best_distance = distance;
                    if (length == max_length) {
                        break;
                    }
                }
                const std::size_t next = prev[candidate % kWindow];
                if (next == kNone || next >= candidate) {
                    break;
                }
                candidate = next;
            }
        }

        if (best_length >= kMinMatch) {
            write_match(writer, static_cast<unsigned>(best_length),
                        static_cast<unsigned>(best_distance));
            for (std::size_t k = 0; k < best_length; ++k) {
                insert(i + k);
            }
            i += best_length;
        } else {
            write_literal_length_symbol(writer, data[i]);
            insert(i);
            ++i;
        }
    }

    write_literal_length_symbol(writer, 256); // end of block
    writer.flush();
}

void write_chunk(Bytes& out, std::string_view type, const Bytes& payload) {
    put_u32_be(out, static_cast<std::uint32_t>(payload.size()));
    const std::size_t type_offset = out.size();
    out.insert(out.end(), type.begin(), type.end());
    out.insert(out.end(), payload.begin(), payload.end());
    put_u32_be(out, crc32(out.data() + type_offset, out.size() - type_offset));
}

std::uint8_t paeth(int a, int b, int c) {
    const int p = a + b - c;
    const int pa = std::abs(p - a);
    const int pb = std::abs(p - b);
    const int pc = std::abs(p - c);
    if (pa <= pb && pa <= pc) {
        return static_cast<std::uint8_t>(a);
    }
    return static_cast<std::uint8_t>(pb <= pc ? b : c);
}

// Applies the PNG filter that minimizes the sum of absolute (signed) residuals for each row,
// which is the heuristic recommended by the PNG specification.
Bytes filter_scanlines(const Image& image) {
    constexpr std::size_t kBpp = 3;
    const std::size_t stride = image.width() * kBpp;
    Bytes filtered;
    filtered.reserve((stride + 1) * image.height());

    Bytes previous(stride, 0);
    Bytes current(stride);
    std::array<Bytes, 5> candidates;
    for (auto& candidate : candidates) {
        candidate.resize(stride);
    }

    for (std::size_t y = 0; y < image.height(); ++y) {
        const Rgb* row = image.row(y);
        for (std::size_t x = 0; x < image.width(); ++x) {
            current[x * kBpp] = row[x].r;
            current[x * kBpp + 1] = row[x].g;
            current[x * kBpp + 2] = row[x].b;
        }

        std::size_t best_filter = 0;
        unsigned long long best_cost = ~0ull;
        for (std::size_t filter = 0; filter < candidates.size(); ++filter) {
            Bytes& out = candidates[filter];
            unsigned long long cost = 0;
            for (std::size_t i = 0; i < stride; ++i) {
                const int a = i >= kBpp ? current[i - kBpp] : 0;
                const int b = previous[i];
                const int c = i >= kBpp ? previous[i - kBpp] : 0;
                int predictor = 0;
                switch (filter) {
                case 1:
                    predictor = a;
                    break;
                case 2:
                    predictor = b;
                    break;
                case 3:
                    predictor = (a + b) / 2;
                    break;
                case 4:
                    predictor = paeth(a, b, c);
                    break;
                default:
                    break;
                }
                const auto value = static_cast<std::uint8_t>(current[i] - predictor);
                out[i] = value;
                cost += static_cast<unsigned>(std::abs(static_cast<std::int8_t>(value)));
            }
            if (cost < best_cost) {
                best_cost = cost;
                best_filter = filter;
            }
        }

        filtered.push_back(static_cast<std::uint8_t>(best_filter));
        filtered.insert(filtered.end(), candidates[best_filter].begin(),
                        candidates[best_filter].end());
        std::swap(previous, current);
    }
    return filtered;
}

bool has_extension(const std::string& path, std::string_view extension) {
    if (path.size() < extension.size()) {
        return false;
    }
    return std::equal(extension.rbegin(), extension.rend(), path.rbegin(), [](char a, char b) {
        return a == static_cast<char>(b >= 'A' && b <= 'Z' ? b - 'A' + 'a' : b);
    });
}

} // namespace

std::uint32_t crc32(const std::uint8_t* data, std::size_t size, std::uint32_t crc) {
    const auto& table = crc_table();
    crc = ~crc;
    for (std::size_t i = 0; i < size; ++i) {
        crc = table[(crc ^ data[i]) & 0xFFu] ^ (crc >> 8);
    }
    return ~crc;
}

std::uint32_t adler32(const std::uint8_t* data, std::size_t size) {
    constexpr std::uint32_t kMod = 65521;
    std::uint32_t a = 1;
    std::uint32_t b = 0;
    for (std::size_t i = 0; i < size; ++i) {
        a = (a + data[i]) % kMod;
        b = (b + a) % kMod;
    }
    return (b << 16) | a;
}

std::vector<std::uint8_t> zlib_compress(const std::vector<std::uint8_t>& data) {
    Bytes out{0x78, 0x9C}; // deflate, 32K window, default compression
    deflate_fixed(data, out);
    put_u32_be(out, adler32(data.data(), data.size()));
    return out;
}

std::vector<std::uint8_t> encode_png(const Image& image) {
    Bytes out{0x89, 'P', 'N', 'G', '\r', '\n', 0x1A, '\n'};

    Bytes header;
    put_u32_be(header, static_cast<std::uint32_t>(image.width()));
    put_u32_be(header, static_cast<std::uint32_t>(image.height()));
    header.push_back(8); // bit depth
    header.push_back(2); // color type: truecolor
    header.push_back(0); // compression: deflate
    header.push_back(0); // filter method: adaptive
    header.push_back(0); // no interlacing
    write_chunk(out, "IHDR", header);
    write_chunk(out, "IDAT", zlib_compress(filter_scanlines(image)));
    write_chunk(out, "IEND", {});
    return out;
}

std::vector<std::uint8_t> encode_ppm(const Image& image) {
    const std::string header =
        "P6\n" + std::to_string(image.width()) + " " + std::to_string(image.height()) + "\n255\n";
    Bytes out(header.begin(), header.end());
    out.reserve(out.size() + image.width() * image.height() * 3);
    for (std::size_t y = 0; y < image.height(); ++y) {
        const Rgb* row = image.row(y);
        for (std::size_t x = 0; x < image.width(); ++x) {
            out.push_back(row[x].r);
            out.push_back(row[x].g);
            out.push_back(row[x].b);
        }
    }
    return out;
}

bool write_image(const std::string& path, const Image& image) {
    const Bytes data = has_extension(path, ".ppm") ? encode_ppm(image) : encode_png(image);
    std::ofstream file(path, std::ios::binary);
    if (!file) {
        return false;
    }
    file.write(reinterpret_cast<const char*>(data.data()),
               static_cast<std::streamsize>(data.size()));
    return static_cast<bool>(file);
}

} // namespace mandelbrot
