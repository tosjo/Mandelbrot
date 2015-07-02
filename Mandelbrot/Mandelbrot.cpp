// mandelbrot.cpp
// compile with: g++ -std=c++11 mandelbrot.cpp -o mandelbrot
// view output with: eog mandelbrot.ppm

#include "stdafx.h"
#include <fstream>
#include <complex> // if you make use of complex number facilities in C++
#include <iostream>
#include <cstdlib>
#include <thread>
#include <mutex>
#include <vector>

using namespace std;

template <class T> struct RGB { T r, g, b; };

template <class T>
struct Matrix
{
	std::vector<T> data;
	size_t rows;
	size_t cols;

	class proxy {
		Matrix &m;
		size_t index_1;
	public:
		proxy(Matrix &m, size_t index_1) : m(m), index_1(index_1) { }

		T &operator[](size_t index) { return m.data[index * m.rows + index_1]; }
	};

	class const_proxy {
		Matrix const &m;
		size_t index_1;
	public:
		const_proxy(Matrix const &m, size_t index_1) : m(m), index_1(index_1) { }

		T const &operator[](size_t index) const { return m.data[index * m.rows + index_1]; }
	};


public:
	Matrix(size_t rows, size_t cols) : data(rows * cols), rows(rows), cols(cols) { }

	proxy operator[](size_t index) { return proxy(*this, index); }
	const_proxy operator[](size_t index) const { return const_proxy(*this, index); }

};

template <class T>
std::ostream &operator<<(std::ostream &out, Matrix<T> const &m) {
	out << "P6" << std::endl << m.cols << " " << m.rows << std::endl << 255 << std::endl;
	for (size_t y = 0; y < m.rows; y++)
	for (size_t x = 0; x < m.cols; x++) {
		T pixel = m[y][x];
		out << pixel.r << pixel.g << pixel.b;
	}
	return out;
}

/*Draw Mandelbrot according to the provided parameters*/
template <class T>
void draw_Mandelbrot(T & image, const unsigned width, const unsigned height, double cxmin, double cxmax, double cymin, double cymax, unsigned int max_iterations) {

#pragma omp parallel for
	for (int ix = 0; ix < width; ++ix)
	for (int iy = 0; iy < height; ++iy)
	{
		std::complex<double> c(cxmin + ix / (width - 1.0)*(cxmax - cxmin), cymin + iy / (height - 1.0)*(cymax - cymin));
		std::complex<double> z = 0;
		unsigned int iterations;

		for (iterations = 0; iterations < max_iterations && std::abs(z) < 2.0; ++iterations)
			z = z*z + c;

		image[iy][ix].r = image[iy][ix].g = image[iy][ix].b = iterations;

	}
}

int main() {
	const unsigned width = 1600;
	const unsigned height = 1600;

	Matrix<RGB<unsigned char>> image(height, width);

	clock_t start = clock();
	draw_Mandelbrot(image, width, height, -2.0, 0.5, -1.0, 1.0, 255);
	clock_t stop = clock();

	std::cout << (double(stop - start) / CLOCKS_PER_SEC) << " seconds\n";

	std::ofstream out("mandelbrot.ppm", std::ios::binary);
	out << image;

	return 0;
}