#include <iostream>
#include <fstream>
#include "Graph.h"
#include "IDLA.h"


int main() {

    std::cout << "Hello, World!" << std::endl;


    std::cout << "starting...." << std::endl;
    auto start = std::chrono::steady_clock::now();

// ----------   RUN IDLA   --------------

    int N = 1024;

    int H = N / 4 ;
    // arguments: (N, H, k)
    // := (cylinder width, height of interface for amount of time, k value for e^(i k x) eigenfunction)
    IDLA graph(N, H, 1);
    graph.initializeGraph();
    // second argument here must be N * H so that expected interface height is H for each sample
    graph.MarkovChain(100, N * H );



    // OPEN FILE

    ofstream myfile;
    myfile.open ("IDLA_DATA/corr9.txt");


    // write correlations to file
    for (auto c : graph.correlations) {
        if (abs(c) < .00001) {
            myfile << 0;
        } else {
            myfile << c;
        }
        myfile << ", ";
    }

    myfile.close();


    ofstream myfile2;
    myfile2.open ("IDLA_DATA/fourier9.txt");

    for (auto f : graph.Fourier) {
        if (std::abs(real(f)) < 0.0001)
            f.real(0);
        if (std::abs(imag(f)) < 0.0001)
            f.imag(0);
        myfile2 << real(f) << " + I * ("<< imag(f) << ")";
        myfile2 << ", ";
    }

    myfile2.close();

    ofstream myfile5;
    myfile5.open ("IDLA_DATA/Dphi9.txt");

    // write end martingale values to file
    for (complex<double> m : graph.Dphivals) {
        if (std::abs(real(m)) < 0.0001)
            m.real(0);
        if (std::abs(imag(m)) < 0.0001)
            m.imag(0);
        myfile5 << real(m) << " + I * ("<< imag(m) << ")";
        myfile5 << ", ";
    }
    myfile5.close();


    ofstream myfile4;
    myfile4.open ("IDLA_DATA/mart9.txt");

    // write end martingale values to file
    for (complex<double> m : graph.Mvals) {
        myfile4 << real(m) << " + I * ("<< imag(m) << ")";
        myfile4 << ", ";
    }
    myfile4.close();

    cout << graph.Fourier.size() << endl;

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

//    int N = 512;
//
//    // arguments: (X, Y, k)
//    // := (cylinder width, cylinder height, k value for e^(i k x) eigenfunction)
//    CylinderGraph graph(N, N, 3);
//    graph.initializeGraph();
//
//    // Here, we will try to sample Mt from the same chain but at far apart times to decorrelate samples
//    graph.MarkovChain(N, N);
//
//
//
//    // OPEN FILE
//
//    ofstream myfile;
//    myfile.open ("EROSION_DATA/M3.txt");



//    // write end martingale values to file
//    for (complex<double> m : graph.Mvals) {
//        myfile << real(m) << " + I * ("<< imag(m) << ")";
//        myfile << ", ";
//    }


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

    // write martingale to file
//    for (complex<double> m : graph.M) {
//        myfile << real(m) << " + I * ("<< imag(m) << ")";
//        myfile << ", ";
//    }

//    myfile.close();


    // graph.printGraph();

    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    typedef std::chrono::duration<float> float_seconds;
    auto secs = std::chrono::duration_cast<float_seconds>(diff);
    std::cout << "Time elapsed: " << secs.count() << " s." << std::endl;






    return 0;


}