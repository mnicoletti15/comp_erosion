#include <iostream>
#include <fstream>
#include "Graph.h"
#include "IDLA.h"


int main() {

    std::cout << "Hello, World!" << std::endl;


    std::cout << "starting...." << std::endl;
    auto start = std::chrono::steady_clock::now();

// ----------   RUN IDLA   --------------
//
//    int N = 1024;
//
//    int H = N / 4 ;
//    // arguments: (N, H, k)
//    // := (cylinder width, height of interface for amount of time, k value for e^(i k x) eigenfunction)
//    IDLA graph(N, H, 1);
//    graph.initializeGraph();
//    // second argument here must be N * H so that expected interface height is H for each sample
//    graph.MarkovChain(100, N * H );
//
//
//
//    // OPEN FILE
//
//    ofstream myfile;
//    myfile.open ("IDLA_DATA/corr9.txt");
//
//
//    // write correlations to file
//    for (auto c : graph.correlations) {
//        if (abs(c) < .00001) {
//            myfile << 0;
//        } else {
//            myfile << c;
//        }
//        myfile << ", ";
//    }
//
//    myfile.close();
//
//
//    ofstream myfile2;
//    myfile2.open ("IDLA_DATA/fourier9.txt");
//
//    for (auto f : graph.Fourier) {
//        if (std::abs(real(f)) < 0.0001)
//            f.real(0);
//        if (std::abs(imag(f)) < 0.0001)
//            f.imag(0);
//        myfile2 << real(f) << " + I * ("<< imag(f) << ")";
//        myfile2 << ", ";
//    }
//
//    myfile2.close();
//
//    ofstream myfile5;
//    myfile5.open ("IDLA_DATA/Dphi9.txt");
//
//    // write end martingale values to file
//    for (complex<double> m : graph.Dphivals) {
//        if (std::abs(real(m)) < 0.0001)
//            m.real(0);
//        if (std::abs(imag(m)) < 0.0001)
//            m.imag(0);
//        myfile5 << real(m) << " + I * ("<< imag(m) << ")";
//        myfile5 << ", ";
//    }
//    myfile5.close();
//
//
//    ofstream myfile4;
//    myfile4.open ("IDLA_DATA/mart9.txt");
//
//    // write end martingale values to file
//    for (complex<double> m : graph.Mvals) {
//        myfile4 << real(m) << " + I * ("<< imag(m) << ")";
//        myfile4 << ", ";
//    }
//    myfile4.close();
//


//    ofstream myfile3;
//    myfile3.open ("IDLA_DATA/mart.txt");
//
//    cout << graph.M.size() << endl;
//
//    for (complex<double> Mt : graph.M) {
//        if (std::abs(real(Mt)) < 0.0001)
//            Mt.real(0);
//        if (std::abs(imag(Mt)) < 0.0001)
//            Mt.imag(0);
//        myfile3 << real(Mt) << " + I * ("<< imag(Mt) << ")";
//        myfile3 << ", ";
//    }
//
//    myfile3.close();




//   ----------------------------------------- End of IDLA --------------------------------------




    //    RUN EROSION

//    int N = 256;
//    int t0 = N;
//    int k = 3;
//    int T0 = t0 * N;
//
//    // arguments: (X, Y, k)
//    // := (cylinder width, cylinder height, k value for e^(i k x) eigenfunction)
//    CylinderGraph graph(N, N, k, T0);
//    graph.initializeGraph();
//
//
//
//    graph.MarkovChain(t0, N);

// Run with this setup to get many samples at same time T = N * N
    int N = 256;
    int k = 3;
    int T0 = N * N;

    // arguments: (X, Y, k)
    // := (cylinder width, cylinder height, k value for e^(i k x) eigenfunction)
    CylinderGraph graph(N, N, k, T0);
    graph.initializeGraph();



    graph.MarkovChain(100, N * N);



    // OPEN FILE
//
    ofstream myfile1;
    myfile1.open ("EROSION_DATA/data3/M.txt");

    // write martingale trajectory values to file
    for (complex<double> m : graph.M) {
        if (std::abs(real(m)) < 0.0001)
            m.real(0);
        myfile1 << real(m);

        if (std::abs(imag(m)) >= 0.0001)
            cout << " + I * ("<< imag(m) << ")";

        myfile1 << ", ";
    }


    myfile1.close();
//
//    ofstream myfile2;
//    myfile2.open ("EROSION_DATA/data2/MQ.txt");
//
//    //     write martingale difference values to file
//
//    for (double m : graph.MQ) {
//        if (std::abs(m) < 0.0001)
//            myfile2 << 0;
//        else {
//            myfile2 << m;
//        }
//
//        myfile2 << ", ";
//    }
//
//    myfile2.close();
//
//
//
    ofstream myfile3;
    myfile3.open ("EROSION_DATA/data3/Dphi.txt");


    // write Dphi values to file
    for (complex<double> m : graph.Dphilist) {
        if (std::abs(real(m)) < 0.0001)
            m.real(0);
        myfile3 << real(m);

        if (std::abs(imag(m)) >= 0.0001)
            cout << " + I * ("<< imag(m) << ")";

        myfile3 << ", ";
    }
    myfile3.close();
//
//
//    ofstream myfile4;
//    myfile4.open ("EROSION_DATA/data2/DphiQ.txt");
//
//    //     write Dphi difference values to file
//
//    for (double m : graph.DphiQ) {
//        if (std::abs(m) < 0.0001) {
//            myfile4 << 0;
//        } else {
//            myfile4 << m;
//        }
//
//        myfile4 << ", ";
//    }
//
//    myfile4.close();
////
////
//    ofstream myfile5;
//    myfile5.open ("EROSION_DATA/data2/Y.txt");
//
//    //     write Y values to file
//
//    for (complex<double> m : graph.Ylist) {
//        if (std::abs(real(m)) < 0.0001)
//            m.real(0);
//        myfile5 << real(m);
//
//        if (std::abs(imag(m)) >= 0.0001)
//            cout << " + I * ("<< imag(m) << ")";
//
//        myfile5 << ", ";
//    }
//
//    myfile5.close();
////
////
//    ofstream myfile6;
//    myfile6.open ("EROSION_DATA/data2/YQ.txt");
//
//    //     write Y quad var values to file
//
//    for (double m : graph.YQ) {
//        if (std::abs(m) < 0.0001) {
//            myfile6 << 0;
//        } else {
//            myfile6 << m;
//        }
//
//        myfile6 << ", ";
//    }
//
//    myfile6.close();


    /*
    // write fluctuations to file
    for (auto x : graph.fluctuations) {
        cout << x << ", ";
    }
    std::cout << std::endl;

    // write correlations to file
    for (auto c : graph.correlations) {
        myfile << c;
        myfile << ", ";
    }
    */

    int L1 = graph.MQ.size();
    cout << "final value for MQ: " << graph.MQ[L1 - 1] << endl;
    cout << "slope: " << graph.MQ[L1 - 1] / T0 * N * N << endl;


    int L2 = graph.DphiQ.size();
    cout << "final value for DphiQ: " << graph.DphiQ[L2 - 1] << endl;
    cout << "slope: " << graph.DphiQ[L2 - 1] / T0 * N * N << endl;


    int L3 = graph.YQ.size();
    cout << "final value Y quad var: " << graph.YQ[L3 - 1] << endl;

//     graph.printGraph();

    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    typedef std::chrono::duration<float> float_seconds;
    auto secs = std::chrono::duration_cast<float_seconds>(diff);
    std::cout << "Time elapsed: " << secs.count() << " s." << std::endl;






    return 0;


}