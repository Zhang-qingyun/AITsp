#include "Simulator.h"


using namespace std;
using namespace szx;


int main() {
    //Simulator::initDefaultEnvironment();

    Simulator sim;
    sim.debug();
    //sim.benchmark(1);
    //sim.parallelBenchmark(1);
    //sim.generateInstance();
   // for (int i = 1; i <= 40; ++i) { sim.convertPmedInstance("Instance/pmed/pmed", i); }
   /* sim.convertTspInstance("d493", 10);
    sim.convertTspInstance("d657", 10);
    sim.convertTspInstance("gr202", 10);
    sim.convertTspInstance("kroA200", 10);
    sim.convertTspInstance("kroB200", 10);
    sim.convertTspInstance("pcb442", 10);
    sim.convertTspInstance("pr226", 10);
    sim.convertTspInstance("pr264", 10); 
    sim.convertTspInstance("pr299", 10); */
   // for (int p = 10; p <= 150; p += 10) { sim.convertTspInstance("u1060", p); }
   // for (int p = 10; p <= 150; p += 10) { sim.convertTspInstance("pcb3038", p); }
   // for (int p = 10; p <= 150; p += 10) { sim.convertTspInstance("u1817", p); }
   // for (int p = 10; p <= 100; p += 10) { sim.convertTspInstance("rl1323", p); }
   // for (int p = 10; p <= 50; p += 10) { sim.convertTspInstance("pcb3038", p); }
   // for (int p = 100; p <= 500; p += 50) { sim.convertTspInstance("pcb3038", p); }

    return 0;
}
