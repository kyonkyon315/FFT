#pragma once
#include <string>
#include <iostream>
#include <vector>
#include <complex>
#include <stdexcept>
#include <fstream>
class ComplexData
{
private:
	double* m_real;
	double* m_imag;
	size_t m_dataSize;
public:
	ComplexData(size_t dataSize);
	~ComplexData();
	void loadtxt(std::string filename);
	void savetxt(std::string filename);
	void setValueZero();
	size_t dataSize() const { return m_dataSize; }
	double* data_real() { return m_real; }
	double* data_imag() { return m_imag; }
	double& operator[](size_t i) { return m_real[i]; }
};


