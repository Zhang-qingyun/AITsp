#include "Simulator.h"


using namespace std;
using namespace zqy;


int main() {
    //Simulator::initDefaultEnvironment();

    Simulator sim;
    while(1)
    sim.debug();
    //sim.benchmark(1);

    return 0;
}
