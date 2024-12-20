#include "ComplexData.h"

ComplexData::ComplexData(size_t dataSize) 
	:m_real(nullptr)
	,m_imag(nullptr)
	,m_dataSize()
{
	if (dataSize == 0) {
		std::cerr << "dataSize must be larger than 0.\n";
		throw std::runtime_error("");
	}
	try {
		m_real = new double[dataSize];
		m_imag = new double[dataSize];
	}
	catch (...) {
		delete[] m_real;
		delete[] m_imag;
		throw;
	}
	m_dataSize = dataSize;
	for (size_t i = 0; i < m_dataSize; i++) {
		m_real[i] = 0.;
		m_imag[i] = 0.;
	}
}

ComplexData::~ComplexData() {
	delete[] m_real;
	delete[] m_imag;
}

void ComplexData::loadtxt(std::string filename) {
	std::ifstream ifs(filename);
	for (size_t i = 0; i < m_dataSize; i++) {
		ifs>>m_real[i];
	}
}

void ComplexData::savetxt(std::string filename) {
	std::ofstream ofs(filename);
	ofs << m_real[0];
	for (size_t i = 1; i < m_dataSize; i++) {
		ofs <<" "<< m_real[i];
	}
	ofs << "\n"<<m_imag[0];
	for (size_t i = 1; i < m_dataSize; i++) {
		ofs <<" "<< m_imag[i];
	}
}

void ComplexData::setValueZero() {
	for (size_t i = 1; i < m_dataSize; i++) {
		m_real[i] = 0.;
		m_imag[i] = 0.;
	}
}


