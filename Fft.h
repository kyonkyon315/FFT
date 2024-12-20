#pragma once
#include "ComplexData.h"
#define U64 unsigned long long
class Fft
{
private:
	bool isIfft;

	size_t m_dataSize;

	U64 i, j, k;
	long long l, m, n;

	U64 maxDepth;

	double* table_r;
	double* table_i;
	double* w_r;
	double* w_i;

	U64 currentDepth ;
	char* childStatus;//'a':childEven childOddã§Ç…ñ¢äÆê¨, 'b':childEvenÇÃÇ›äÆê¨, 'c':ã§Ç…äÆê¨
	double** iputPoint_r;
	double** iputPoint_i;
	double** oputPoint_r;
	double** oputPoint_i;


	double* currentIputPoint_r, * currentOputPoint_r, * currentIputChild_r, * currentOputChild_r;
	double* currentIputPoint_i, * currentOputPoint_i, * currentIputChild_i, * currentOputChild_i;
	double* currentW_r;
	double* currentW_i;
public:
	Fft();
	~Fft();

	void setDataSize(U64 dataSize);
	void fft(ComplexData& input, ComplexData& output, bool idft = false);

};

