#pragma once
#include <vector>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include "ComplexData.h"
class Dft
{
private:
	size_t m_dataSize;
	double* w_real;
	double* w_imag;
	double angle;

	bool isIdft;
public:
	Dft();
	Dft(const size_t dataSize,bool idft=false);
	~Dft();

	
	void dft(ComplexData& input,ComplexData& output,const size_t dataSize,bool idft=false);
};


