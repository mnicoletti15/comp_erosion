#include <iostream>
#include <fstream>
#include "Graph.h"


int main() {

    std::cout << "Hello, World!" << std::endl;


    std::cout << "starting...." << std::endl;
    auto start = std::chrono::steady_clock::now();

    CylinderGraph graph(64, 64);
    graph.initializeGraph();
    graph.MarkovChain(50, 400);


    ofstream myfile;
    myfile.open ("example.txt");


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
    for (complex<double> m : graph.M) {
        myfile << real(m) << " + I * ("<< imag(m) << ")";
        myfile << ", ";
    }

    myfile.close();

    graph.printGraph();

    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    typedef std::chrono::duration<float> float_seconds;
    auto secs = std::chrono::duration_cast<float_seconds>(diff);
    std::cout << "Time elapsed: " << secs.count() << " s." << std::endl;

    return 0;


}