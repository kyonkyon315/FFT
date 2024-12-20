#include <iostream>
#include "Dft.h"
#include "ComplexData.h"
#include <fstream>
#include "Timer.h"
#include "Fft.h"

int main()
{
    Fft fft;
    ComplexData ew(36000);
    ComplexData ns(36000);
    ComplexData ud(36000);

    //txtファイルの読み取り
    ew.loadtxt("ew.txt");
    ns.loadtxt("ns.txt");
    ud.loadtxt("ud.txt");

    ComplexData absData(36000);

    //振幅の絶対値を計算
    for (int i = 0; i < 36000; i++) {
        absData[i] = sqrt(ew[i] * ew[i] + ns[i] * ns[i] + ud[i] * ud[i]);
    }

    Timer timer;

    //フーリエ変換の結果を入れるデータ
    ComplexData Y(36000);


    for (int i = 0; i <= 15; i++) {
        
        fft.setDataSize(1ULL << i);//dataSizeを2^iに設定
        Y.setValueZero();

        timer.start();//計測開始

        fft.fft(absData, Y);//高速フーリエ変換の実行

        timer.stop();//計測終了

        std::cout << "N: 2^" << std::setw(3) << i 
                  << " -> " 
                  << timer << "\n";
        
        Y.savetxt("output" + std::to_string(i) + ".txt");//計算結果を保存
    }

    return 0;
}


