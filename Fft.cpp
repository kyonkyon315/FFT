#include "Fft.h"
#define U64 unsigned long long
#define M_PI 3.14159265358979323846

Fft::~Fft() {
	delete[] childStatus;
	delete[] iputPoint_r;
	delete[] iputPoint_i;
	delete[] oputPoint_r;
	delete[] oputPoint_i;
	delete[] table_r;
	delete[] table_i;
	delete[] w_r;
	delete[] w_i;
}

Fft::Fft()
	: childStatus(nullptr), iputPoint_r(nullptr), iputPoint_i(nullptr)
	, oputPoint_r(nullptr), oputPoint_i(nullptr), m_dataSize()
	, table_r(nullptr), table_i(nullptr), w_r(nullptr), w_i(nullptr)
	, currentIputChild_r(), currentIputChild_i(), currentOputChild_r(), currentOputChild_i()
	, currentIputPoint_r(), currentIputPoint_i(), currentOputPoint_r(), currentOputPoint_i()
	, currentDepth(), currentW_r(), currentW_i()
	, i(), j(), k()
	, l(), m(), n()
	, maxDepth(), isIfft()
{
	try {
		childStatus = new char[40];//'a':childEven childOdd共に未完成, 'b':childEvenのみ完成, 'c':共に完成
		iputPoint_r = new double* [40];
		iputPoint_i = new double* [40];
		oputPoint_r = new double* [40];
		oputPoint_i = new double* [40];
	}
	catch (...) {
		delete[]childStatus;
		delete[]iputPoint_r;
		delete[]iputPoint_i;
		delete[]oputPoint_r;
		delete[]oputPoint_i;
		childStatus = nullptr;
		iputPoint_r = nullptr;
		iputPoint_i = nullptr;
		oputPoint_r = nullptr;
		oputPoint_i = nullptr;
		throw;
	}
}

void Fft::setDataSize(U64 dataSize) {
	if (m_dataSize == dataSize)return;
	delete[] table_r;
	delete[] table_i;
	delete[] w_r;
	delete[] w_i;
	table_r = nullptr;
	table_i = nullptr;
	w_r = nullptr;
	w_i = nullptr;
	for(i = 63ULL; i >= 0ULL; i--) {
		if (dataSize & (1ULL << i)) break;
	}
	if (dataSize != (1ULL << i)) {
		std::cout << " Alart: dataSize must be power of two.\n"
			<< " dataSize ( = " << dataSize << " ) is changed to 2^" << i << " ( = ";
		dataSize = (1ULL << i);
		std::cout << dataSize << " ).\n";
	}
	maxDepth = i;// pow(2, maxDepth) = dataSize
	if (maxDepth == 0)return;
	double angleUnit = M_PI / (double)(dataSize >> 1);//error
	double angle = 0.;

	try {
		table_r = new double[dataSize << 1];//データ長の二倍確保すりゅ。
		table_i = new double[dataSize << 1];//データ長の二倍確保すりゅ。
		w_r = new double[dataSize >> 1];
		w_i = new double[dataSize >> 1];
	}
	catch (...) {
		delete[] table_r;
		delete[] table_i;
		delete[] w_r;
		delete[] w_i;
		std::cerr << " Alart: at Fft::setDataSize(dataSize) in " __FILE__ "\n"
			" problem allocating memory.\n";
		throw;
	}
	m = dataSize >> 1;//error
	for (i = 0ULL; i < m; i++) {
		w_r[i] = cos(angle);
		w_i[i] = sin(angle);
		angle += angleUnit;
	}//w[k]=exp{ j * pi k/(dataSize/2) }  jは虚数単位


	/**********************************************************************************
	* キャッシュヒット率向上を期待して、メモリパターンを以下のようにする。            *
	* x0,x1,...xn( = Yn ),Yn-1,...y1,y0                                               *
	* ただし、xkはバタフライ演算前半のk回目に使うメモリ                               *
	* ただし、Ykはバタフライ演算後半の(n-k)回目に使うメモリ                           *
	***********************************************************************************/
	for (U64 depth = 1ULL; depth <= maxDepth; depth++) {
		iputPoint_r[depth] = table_r + dataSize - (dataSize >> (depth - 1ULL));
		iputPoint_i[depth] = table_i + dataSize - (dataSize >> (depth - 1ULL));
		oputPoint_r[depth] = table_r + dataSize + (dataSize >> (depth       )) - 3LL;
		oputPoint_i[depth] = table_i + dataSize + (dataSize >> (depth       )) - 3LL;
	}//iputPoint[maxDepth]==oputPoint[maxDepth]が成り立つようにできている。
	m_dataSize = dataSize;
}

void Fft::fft(ComplexData& input, ComplexData& output, bool idft) {
	if (input.dataSize() < m_dataSize || input.dataSize() < m_dataSize) {
		std::cerr << "Alart: at Fft::fft in" __FILE__ "\n"
			"input.dataSize() and output.dataSize must be larger than this->m_dataSize \n"
			"in order to resize this->m_dataSize, call Fft::setDataSize()\n"<<std::flush;
		throw std::runtime_error("Fft::fft()");
	}
	//深さ0のバタフライ演算のメモリは引数ComplexData内の配列にする。
	iputPoint_r[0] = input.data_real();
	iputPoint_i[0] = input.data_imag();
	oputPoint_r[0] = output.data_real();
	oputPoint_i[0] = output.data_imag();
	if (maxDepth == 0) {
		oputPoint_r[0][0] = iputPoint_r[0][0];
		oputPoint_i[0][0] = iputPoint_i[0][0];
		return;
	}
	/*
	 バタフライ演算の挙動を
	 'a':childEven childOdd共に未完成, 
	 'b':childEvenのみ完成, 
	 'c':共に完成
	 の3状態で管理する*/

	//始めの状態は'a'
	childStatus[0] = 'a';

    //始めは深さ0
	currentDepth = 0;

	for (;;) {
		//even,odd共にfftできていない場合
		if (childStatus[currentDepth]=='a') {
			if (currentDepth == maxDepth) {
				currentDepth--;
				continue;
			}
			currentIputPoint_r = iputPoint_r[currentDepth];
			currentIputPoint_i = iputPoint_i[currentDepth];
			currentIputChild_r = iputPoint_r[currentDepth + 1];
			currentIputChild_i = iputPoint_i[currentDepth + 1];

			j = m_dataSize >> (currentDepth + 1ULL);//=childの長さ
			for (i = 0ULL; i < j; i++) {
				currentIputChild_r[i] = currentIputPoint_r[i << 1];
				currentIputChild_i[i] = currentIputPoint_i[i << 1];
			}
			childStatus[currentDepth] = 'b';
			childStatus[currentDepth + 1] = 'a';
			currentDepth++;
		}

		//even,のみfftできている場合
		else if (childStatus[currentDepth] == 'b') {
			currentIputPoint_r = iputPoint_r[currentDepth];
			currentIputPoint_i = iputPoint_i[currentDepth];
			currentOputPoint_r = oputPoint_r[currentDepth];
			currentOputPoint_i = oputPoint_i[currentDepth];
			currentIputChild_r = iputPoint_r[currentDepth + 1];
			currentIputChild_i = iputPoint_i[currentDepth + 1];
			currentOputChild_r = oputPoint_r[currentDepth + 1];
			currentOputChild_i = oputPoint_i[currentDepth + 1];

			j = m_dataSize >> (currentDepth + 1ULL);//childの長さ
			for (i = 0ULL; i < j; i++) {
				currentOputPoint_r[i] = currentOputChild_r[i];
				currentOputPoint_i[i] = currentOputChild_i[i];
				currentOputPoint_r[j+i] = currentOputChild_r[i];
				currentOputPoint_i[j+i] = currentOputChild_i[i];
			}
			for (i = 0ULL; i < j; i++) {
				currentIputChild_r[i] = currentIputPoint_r[(i << 1) + 1ULL];
				currentIputChild_i[i] = currentIputPoint_i[(i << 1) + 1ULL];
			}

			childStatus[currentDepth] = 'c';
			childStatus[currentDepth + 1] = 'a';
			currentDepth++;
		}

		//even,oddともにfftできている場合
		else {
			currentOputPoint_r = oputPoint_r[currentDepth];
			currentOputPoint_i = oputPoint_i[currentDepth];
			currentOputChild_r = oputPoint_r[currentDepth + 1];
			currentOputChild_i = oputPoint_i[currentDepth + 1];

			currentW_r = w_r;
			currentW_i = w_i;

			j = m_dataSize >> (currentDepth + 1ULL);//childの長さ
			k = (1ULL << currentDepth);
			
			for (i = 0ULL; i < j; i++) {
				currentOputPoint_r[i] += (*currentW_r * currentOputChild_r[i] - (*currentW_i) * currentOputChild_i[i]);
				currentOputPoint_i[i] += (*currentW_r * currentOputChild_i[i] + (*currentW_i) * currentOputChild_r[i]);
				currentOputPoint_r[j+i] -= (*currentW_r * currentOputChild_r[i] - (*currentW_i) * currentOputChild_i[i]);
				currentOputPoint_i[j+i] -= (*currentW_r * currentOputChild_i[i] + (*currentW_i) * currentOputChild_r[i]);
				currentW_r += k;
				currentW_i += k;
					/*
					*currentW= w[i*k]
						 = exp{j M_PI * i * (1ULL << currentDepth)/ (dataSize >> 1)
						 = exp{j M_PI * i / (dataSize>>(currentDepth+1))
						 = exp{j i * M_PI / childの長さ
						 = exp{j 2.* i * M_PI / N}
							  jは虚数単位*/
			}
			if (currentDepth == 0)break;
			currentDepth--;
		}
	}
}
#undef M_PI




