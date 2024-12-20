#include "Dft.h"
#define M_PI 3.14159265358979323846
Dft::Dft()
	:m_dataSize()
	, angle()
	, w_real(nullptr)
	, w_imag(nullptr)
	, isIdft(true)
{
}

Dft::Dft(const size_t dataSize, bool idft) 
	:m_dataSize()
	, angle()
	, w_real(nullptr)
	, w_imag(nullptr)
	, isIdft(true)
{
	if (dataSize == 0U) {
		std::cout << " Alart at Dft::Dft(size_t dataSize, bool idft) in " __FILE__ "\n"
		             " dataSize must be larger than 0\n";
		throw std::runtime_error("Dft::Dft(size_t dataSize, bool idft)");
	}
	try {
		w_real = new double[dataSize];
		w_imag = new double[dataSize];
	}
	catch (...) {
		this->~Dft();
		throw;
	}
	m_dataSize = dataSize;
	angle = 2. * M_PI / static_cast<double>(dataSize)*(idft? -1.:1.);
	isIdft = idft;
	for (size_t i = 0; i < dataSize; ++i) {
		w_real[i] = cos(angle * (double)i);
		w_imag[i] = sin(angle * (double)i);
	}
}

Dft::~Dft() {
	delete[] w_real;
	delete[] w_imag;
	w_real = nullptr;
	w_imag = nullptr;
	m_dataSize = 0U;
}


void Dft::dft(ComplexData& input, ComplexData& output,const size_t dataSize, bool idft) {

	if (dataSize == 0 || dataSize > input.dataSize() || input.dataSize() != output.dataSize()) {
		std::cerr << "Dft::dft()\n";
		throw std::runtime_error("Dft::dft()");
	}
	if (dataSize != m_dataSize) {
		this->~Dft();
		try {
			w_real = new double[dataSize];
			w_imag = new double[dataSize];
		}
		catch (...) {
			this->~Dft();
		}
		m_dataSize = dataSize;
		angle = 2. * M_PI / static_cast<double>(dataSize) * (idft ? -1. : 1.);
		isIdft = idft;
		for (size_t i = 0; i < dataSize; ++i) {
			w_real[i] = cos(angle * (double)i);
			w_imag[i] = sin(angle * (double)i);
		}
	}
	else if (idft != isIdft) {
		for (size_t i = 0; i < m_dataSize; ++i) {
			w_real[i] = -w_real[i];
			w_imag[i] = -w_imag[i];
		}
	}
	double* real_in = input.data_real();
	double* imag_in = input.data_imag();
	double* real_out = output.data_real();
	double* imag_out = output.data_imag();
	size_t nTimesMdivN;
	for (size_t m = 0; m < m_dataSize; m++) {
		nTimesMdivN = 0;
		real_out[m] = 0.;
		imag_out[m] = 0.;
		for (size_t n = 0; n < m_dataSize; n++) {
			//y+=w[nTimesMdivN]*x[n]
			real_out[m] += (w_real[nTimesMdivN] * real_in[n] - w_imag[nTimesMdivN] * imag_in[n]);
			//real_out[m] += (w_real[(m*n)%m_dataSize] * real_in[n] - w_imag[(m * n) % m_dataSize] * imag_in[n]);
			imag_out[m] += (w_real[nTimesMdivN] * imag_in[n] + w_imag[nTimesMdivN] * real_in[n]);
			//imag_out[m] += (w_real[(m*n)%m_dataSize] * imag_in[n] + w_imag[(m * n) % m_dataSize] * real_in[n]);
			nTimesMdivN += m;
			if (nTimesMdivN >= m_dataSize)nTimesMdivN -= m_dataSize;
		}
	}
}

#undef M_PI


