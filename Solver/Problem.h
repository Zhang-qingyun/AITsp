////////////////////////////////
/// usage : 1.	data that identifies a guillotine cut problem and its solution.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef SMART_SZX_P_CENTER_PROBLEM_H
#define SMART_SZX_P_CENTER_PROBLEM_H


#include "Config.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include "Common.h"
#include "PbReader.h"
#include "TSP.pb.h"


namespace zqy {

class Problem {
    #pragma region Type
public:
    struct Input : public pb::TSP::Input {
        bool load(const String &path) { return pb::load(path, *this); }
    };

    struct Output : public pb::TSP::Output {
         ID coverRadius = 0;
    };
    #pragma endregion Type

    #pragma region Constant
public:
    enum {
        MaxDistance = (1 << 28),
        MaxNodeNum = 5000,
        MaxCenterNum = 5000,

        InvalidId = -1,
    };

    static constexpr double TopologicalGraphObjScale = 1;
    static constexpr double GeometricalGraphObjScale = 100;
    #pragma endregion Constant

    #pragma region Constructor
public:
    #pragma endregion Constructor

    #pragma region Method
public:
    static bool isTopologicalGraph(const pb::TSP::Input &input) {
        return input.graph().nodes().empty();
    }
    #pragma endregion Method

    #pragma region Field
public:
    #pragma endregion Field
}; // Problem

}


#endif // SMART_ZQY_P_CENTER_PROBLEM_H
