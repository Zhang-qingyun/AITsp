#include "Solver.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <mutex>

#include <cmath>


using namespace std;


namespace szx {

#pragma region Solver::Cli
int Solver::Cli::run(int argc, char * argv[]) {
    Log(LogSwitch::Szx::Cli) << "parse command line arguments." << endl;
    Set<String> switchSet;
    Map<String, char*> optionMap({ // use string as key to compare string contents instead of pointers.
        { InstancePathOption(), nullptr },
        { SolutionPathOption(), nullptr },
        { RandSeedOption(), nullptr },
        { TimeoutOption(), nullptr },
        { MaxIterOption(), nullptr },
        { JobNumOption(), nullptr },
        { RunIdOption(), nullptr },
        { EnvironmentPathOption(), nullptr },
        { ConfigPathOption(), nullptr },
        { LogPathOption(), nullptr }
        });

    for (int i = 1; i < argc; ++i) { // skip executable name.
        auto mapIter = optionMap.find(argv[i]);
        if (mapIter != optionMap.end()) { // option argument.
            mapIter->second = argv[++i];
        } else { // switch argument.
            switchSet.insert(argv[i]);
        }
    }

    Log(LogSwitch::Szx::Cli) << "execute commands." << endl;
    if (switchSet.find(HelpSwitch()) != switchSet.end()) {
        cout << HelpInfo() << endl;
    }

    if (switchSet.find(AuthorNameSwitch()) != switchSet.end()) {
        cout << AuthorName() << endl;
    }

    Solver::Environment env;
    env.load(optionMap);
    if (env.instPath.empty() || env.slnPath.empty()) { return -1; }

    Solver::Configuration cfg;
    cfg.load(env.cfgPath);

    Log(LogSwitch::Szx::Input) << "load instance " << env.instPath << " (seed=" << env.randSeed << ")." << endl;
    Problem::Input input;
    if (!input.load(env.instPath)) { return -1; }

    Solver solver(input, env, cfg);
    solver.solve();

    pb::Submission submission;
    submission.set_thread(to_string(env.jobNum));
    submission.set_instance(env.friendlyInstName());
    submission.set_duration(to_string(solver.timer.elapsedSeconds()) + "s");

    solver.output.save(env.slnPath, submission);
    #if SZX_DEBUG
    solver.output.save(env.solutionPathWithTime(), submission);
    solver.record();
    #endif // SZX_DEBUG

    return 0;
}
#pragma endregion Solver::Cli

#pragma region Solver::Environment
void Solver::Environment::load(const Map<String, char*> &optionMap) {
    char *str;

    str = optionMap.at(Cli::EnvironmentPathOption());
    if (str != nullptr) { loadWithoutCalibrate(str); }

    str = optionMap.at(Cli::InstancePathOption());
    if (str != nullptr) { instPath = str; }

    str = optionMap.at(Cli::SolutionPathOption());
    if (str != nullptr) { slnPath = str; }

    str = optionMap.at(Cli::RandSeedOption());
    if (str != nullptr) { randSeed = atoi(str); }

    str = optionMap.at(Cli::TimeoutOption());
    if (str != nullptr) { msTimeout = static_cast<Duration>(atof(str) * Timer::MillisecondsPerSecond); }

    str = optionMap.at(Cli::MaxIterOption());
    if (str != nullptr) { maxIter = atoi(str); }

    str = optionMap.at(Cli::JobNumOption());
    if (str != nullptr) { jobNum = atoi(str); }

    str = optionMap.at(Cli::RunIdOption());
    if (str != nullptr) { rid = str; }

    str = optionMap.at(Cli::ConfigPathOption());
    if (str != nullptr) { cfgPath = str; }

    str = optionMap.at(Cli::LogPathOption());
    if (str != nullptr) { logPath = str; }

    calibrate();
}

void Solver::Environment::load(const String &filePath) {
    loadWithoutCalibrate(filePath);
    calibrate();
}

void Solver::Environment::loadWithoutCalibrate(const String &filePath) {
    // EXTEND[szx][8]: load environment from file.
    // EXTEND[szx][8]: check file existence first.
}

void Solver::Environment::save(const String &filePath) const {
    // EXTEND[szx][8]: save environment to file.
}
void Solver::Environment::calibrate() {
    // adjust thread number.
    int threadNum = thread::hardware_concurrency();
    if ((jobNum <= 0) || (jobNum > threadNum)) { jobNum = threadNum; }

    // adjust timeout.
    msTimeout -= Environment::SaveSolutionTimeInMillisecond;
}
#pragma endregion Solver::Environment

#pragma region Solver::Configuration
void Solver::Configuration::load(const String &filePath) {
    // EXTEND[szx][5]: load configuration from file.
    // EXTEND[szx][8]: check file existence first.
}

void Solver::Configuration::save(const String &filePath) const {
    // EXTEND[szx][5]: save configuration to file.
}
#pragma endregion Solver::Configuration

#pragma region Solver
bool Solver::solve() {
    init();

    int workerNum = (max)(1, env.jobNum / cfg.threadNumPerWorker);
    cfg.threadNumPerWorker = env.jobNum / workerNum;
    List<Solution> solutions(workerNum, Solution(this));
    List<bool> success(workerNum);

    Log(LogSwitch::Szx::Framework) << "launch " << workerNum << " workers." << endl;
    List<thread> threadList;
    threadList.reserve(workerNum);
    for (int i = 0; i < workerNum; ++i) {
        // TODO[szx][2]: as *this is captured by ref, the solver should support concurrency itself, i.e., data members should be read-only or independent for each worker.
        // OPTIMIZE[szx][3]: add a list to specify a series of algorithm to be used by each threads in sequence.
        threadList.emplace_back([&, i]() { success[i] = optimize(solutions[i], i); });
    }
    for (int i = 0; i < workerNum; ++i) { threadList.at(i).join(); }

    Log(LogSwitch::Szx::Framework) << "collect best result among all workers." << endl;
    int bestIndex = -1;
    Length bestValue = Problem::MaxDistance;
    for (int i = 0; i < workerNum; ++i) {
        if (!success[i]) { continue; }
        Log(LogSwitch::Szx::Framework) << "worker " << i << " got " << solutions[i].coverRadius << endl;
        if (solutions[i].coverRadius >= bestValue) { continue; }
        bestIndex = i;
        bestValue = solutions[i].coverRadius;
    }

    env.rid = to_string(bestIndex);
    if (bestIndex < 0) { return false; }
    output = solutions[bestIndex];
    return true;
}

void Solver::record() const {
    #if SZX_DEBUG
    int generation = 0;

    ostringstream log;

    System::MemoryUsage mu = System::peakMemoryUsage();

    Length obj = output.coverRadius;
    Length checkerObj = -1;
    bool feasible = check(checkerObj);

    // record basic information.
    log << env.friendlyLocalTime() << ","
        << env.rid << ","
        << env.instPath << ","
        << feasible << "," << (obj - checkerObj) << ",";
    if (Problem::isTopologicalGraph(input)) {
        log << obj << ",";
    } else {
        auto oldPrecision = log.precision();
        log.precision(2);
        log << fixed << setprecision(2) << (obj / aux.objScale) << ",";
        log.precision(oldPrecision);
    }
    log << timer.elapsedSeconds() << ","
        << mu.physicalMemory << "," << mu.virtualMemory << ","
        << env.randSeed << ","
        << cfg.toBriefStr() << ","
        << generation << "," << iteration << ",";

    // record solution vector.
    // EXTEND[szx][2]: save solution in log.
    log << endl;

    // append all text atomically.
    static mutex logFileMutex;
    lock_guard<mutex> logFileGuard(logFileMutex);

    ofstream logFile(env.logPath, ios::app);
    logFile.seekp(0, ios::end);
    if (logFile.tellp() <= 0) {
        logFile << "Time,ID,Instance,Feasible,ObjMatch,Distance,Duration,PhysMem,VirtMem,RandSeed,Config,Generation,Iteration,Solution" << endl;
    }
    logFile << log.str();
    logFile.close();
    #endif // SZX_DEBUG
}

bool Solver::check(Length &checkerObj) const {
    #if SZX_DEBUG
    enum CheckerFlag {
        IoError = 0x0,
        FormatError = 0x1,
        TooManyCentersError = 0x2
    };

    checkerObj = System::exec("Checker.exe " + env.instPath + " " + env.solutionPathWithTime());
    if (checkerObj > 0) { return true; }
    checkerObj = ~checkerObj;
    if (checkerObj == CheckerFlag::IoError) { Log(LogSwitch::Checker) << "IoError." << endl; }
    if (checkerObj & CheckerFlag::FormatError) { Log(LogSwitch::Checker) << "FormatError." << endl; }
    if (checkerObj & CheckerFlag::TooManyCentersError) { Log(LogSwitch::Checker) << "TooManyCentersError." << endl; }
    return false;
    #else
    checkerObj = 0;
    return true;
    #endif // SZX_DEBUG
}

void Solver::init() {
    ID nodeNum = input.graph().nodenum();

    aux.adjMat.init(nodeNum, nodeNum);
    fill(aux.adjMat.begin(), aux.adjMat.end(), Problem::MaxDistance);
    for (ID n = 0; n < nodeNum; ++n) { aux.adjMat.at(n, n) = 0; }

    if (Problem::isTopologicalGraph(input)) {
        aux.objScale = Problem::TopologicalGraphObjScale;
        for (auto e = input.graph().edges().begin(); e != input.graph().edges().end(); ++e) {
            // only record the last appearance of each edge.
            aux.adjMat.at(e->source(), e->target()) = e->length();
            aux.adjMat.at(e->target(), e->source()) = e->length();
        }

        Timer timer(30s);
        constexpr bool IsUndirectedGraph = true;
        IsUndirectedGraph
            ? Floyd::findAllPairsPaths_symmetric(aux.adjMat)
            : Floyd::findAllPairsPaths_asymmetric(aux.adjMat);
        Log(LogSwitch::Preprocess) << "Floyd takes " << timer.elapsedSeconds() << " seconds." << endl;
    } else { // geometrical graph.
        aux.objScale = Problem::GeometricalGraphObjScale;
        for (ID n = 0; n < nodeNum; ++n) {
            double nx = input.graph().nodes(n).x();
            double ny = input.graph().nodes(n).y();
            for (ID m = 0; m < nodeNum; ++m) {
                if (n == m) { continue; }
                aux.adjMat.at(n, m) = static_cast<Length>(aux.objScale * hypot(
                    nx - input.graph().nodes(m).x(), ny - input.graph().nodes(m).y()));
            }
        }
    }

    aux.coverRadii.init(nodeNum);
    fill(aux.coverRadii.begin(), aux.coverRadii.end(), Problem::MaxDistance);
}




bool Solver::init_solution(Solution &sln,ID nodeNum) {//贪心构造
    srand(rand() % MAX);
    int mindictance = MAX;
    vector <int>flag(nodeNum, 0);

    Solution_path.push_back(0);
    flag[Solution_path[0]] = 1;//flag[i]== 1:node i in the Solution_path,visited
    int tempnode1 = Solution_path[0];
    int tempnode2;
    for (int j = 1; j != nodeNum; j++) {
        mindictance = MAX;
        tempnode1 = Solution_path[j - 1];
        for (int i = 1; i != nodeNum; i++) {
            if (aux.adjMat[tempnode1][i] < mindictance && flag[i] == 0) {
                mindictance = aux.adjMat[tempnode1][i];
                tempnode2 = i;
            } else if (aux.adjMat[tempnode1][i] == mindictance && rand() % 2 == 0 && flag[i] == 0) {
                mindictance = aux.adjMat[tempnode1][i];
                tempnode2 = i;
            }
        }
        Solution_path.push_back(tempnode2);
        inpath[Solution_path[j - 1]][Solution_path[j]] = inpath[Solution_path[j]][Solution_path[j - 1]] = 1;
        flag[tempnode2] = 1;
        local_solution += mindictance;
    }
    Solution_path.push_back(0);
    local_solution += aux.adjMat[Solution_path[Solution_path.size() - 2]][0];
    inpath[Solution_path[Solution_path.size() - 2]][Solution_path[Solution_path.size() - 1]] = inpath[Solution_path[Solution_path.size() - 1]][Solution_path[Solution_path.size() - 2]] = 1;
    return 1;
}

bool Solver::optimize(Solution &sln, ID workerId) {
    Log(LogSwitch::Szx::Framework) << "worker " << workerId << " starts." << endl;
    ID nodeNum = input.graph().nodenum();
    ID centerNum = input.centernum();
    int pc=0;
    inpath = vector<vector<int>>(nodeNum, vector<int>(nodeNum, 0));//标记路径中相连的两边，方便找到比相连的交换边
    bool status = true;
    TabuTenure = vector<vector<int>>(nodeNum, vector<int>(nodeNum, 0));
    pair Pair = { -1 };
    bool tempflag = init_solution(sln,nodeNum);//初始解
    best_solution = local_solution;
    if (tempflag == 0) {
        cout << "error!!";
        return 0;
    }
    cout << "init_best_solution" << best_solution << endl;
    /*for (int i = 0; i != nodeNum + 1; i++)
        cout << Solution_path[i] << "\t";
    cout << endl << endl;*/
    /*for (int i = 0; i != nodeNum; i++) {
        for (int j = 0; j != nodeNum; j++) {
            cout << aux.adjMat[i][j] << "\t";
        }
        cout << endl;
    }*/
    clock_t start_time = clock();
    iter = 1;
    while (iter < MAX&& best_solution>50778&& !timer.isTimeOut()) {
        iter++;
        tempflag = find_pair(sln, Pair,nodeNum);
        if (tempflag == 0)
            break;
        change_pair(sln, Pair,nodeNum);//更新交换对    
       /* for (int i = 0; i != nodeNum + 1; i++)
            cout << Solution_path[i] << "\t";
        cout << endl;*/
       /* if (best_solution > local_solution) {
            pc++;
        }*/
        if (best_solution >local_solution) {
            best_solution = local_solution;

        }
        cout << best_solution << "\t";
       /* for (int i = 0; i != nodeNum + 1; i++)
            cout << Solution_path[i] << "\t";
        cout << endl << endl;*/
    }
    for (int i = 0; i != nodeNum + 1; i++)
        cout << Solution_path[i] << "\t";
    cout << endl << "the init best solution:" << best_solution << endl;
    clock_t end_time = clock();
    check(sln,nodeNum);
    /*if (ScInfo.Sc < best_solution)
        best_solution = ScInfo.Sc;*/
   // sln.maxLength = best_solution;
    cout << "the iter: " << iter << endl;
    //cout << "the most solution: " << best_solution << endl;
   // cout << "the true best solution:" << object << endl;
    cout << "the time is:" << (end_time - start_time)*1.0 / CLOCKS_PER_SEC << "s" << endl;
    Log(LogSwitch::Szx::Framework) << "worker " << workerId << " ends." << endl;
    return status;
}
/*
1.初始解：随机的找到一个合法回路（也可以贪心构造）
2.领域动作：交换两个节点，2-opt
3.禁忌动作：禁忌两个节点之间的交换

2-opt+tabu实现
a->b->c->d->e->f
交换不相交的两个节点;b,e
新的边集序列：
a->（e->d->c->b）->f
代价计算：只是a->e对应a->b，b->f对应e->f的代价发生了变化之间的
*/
//路径中node1总是在node2的前面
int Solver::find_pair(Solution &sln, pair &Pair,ID nodeNum) {//确定交换对<nodei，nodej>
    vector <int> id;
    pair tabu_pair;
    tabu_pair.delt = MAX;
    pair no_tabu_pair;
    int delt;
    no_tabu_pair.delt = MAX;
    int no_tabu_samenumber = 1;
    int tabu_samenumber = 1;
    int m = 0;
    for (int i = 1; i != nodeNum - 1; i++) {//起始节点0是不需要变的
        for (int j = i + 1; j != nodeNum; j++) {
            if (inpath[i][j] == 0) {
                //find exchange nodes i,j
                int locali = -1, localj = -1;
                for (int local = 1; local != nodeNum; local++) {
                    if (Solution_path[local] == i)
                        locali = local;

                    if (Solution_path[local] == j)
                        localj = local;
                }
                if (locali < localj) {
                    delt = (aux.adjMat[Solution_path[locali - 1]][j] + aux.adjMat[i][Solution_path[localj + 1]]) - (aux.adjMat[Solution_path[locali - 1]][i] + aux.adjMat[j][Solution_path[localj + 1]]);
                    //计算交换两条边后边的代价变化了多少
                } else {
                    delt = (aux.adjMat[Solution_path[localj - 1]][i] + aux.adjMat[j][Solution_path[locali + 1]]) - (aux.adjMat[Solution_path[localj - 1]][j] + aux.adjMat[i][Solution_path[locali + 1]]);
                }
                if (TabuTenure[i][j] > iter) {//update tabu 
                    if (delt < tabu_pair.delt) {
                        //cout << "nodei,nodej:" << i << "," << j << "," << delt << endl;
                        //cout << "locali,localj:" << locali << "," << localj << endl;
                        if (locali < localj) {
                            tabu_pair.node1 = i;
                            tabu_pair.local1 = locali;
                            tabu_pair.node2 = j;
                            tabu_pair.local2 = localj;
                        } else {
                            tabu_pair.node1 = j;
                            tabu_pair.local1 = localj;
                            tabu_pair.node2 = i;
                            tabu_pair.local2 = locali;
                        }
                        tabu_samenumber = 1;
                        tabu_pair.delt = delt;
                    } else if (delt == tabu_pair.delt) {
                        tabu_samenumber++;
                        if (rand() % tabu_samenumber == 0) {
                            if (locali < localj) {
                                tabu_pair.node1 = i;
                                tabu_pair.local1 = locali;
                                tabu_pair.node2 = j;
                                tabu_pair.local2 = localj;
                            } else {
                                tabu_pair.node1 = j;
                                tabu_pair.local1 = localj;
                                tabu_pair.node2 = i;
                                tabu_pair.local2 = locali;
                            }
                        }
                    }
                } else {//update no_tabu
                    if (delt < no_tabu_pair.delt) {
                       // cout << "nodei,nodej:" <<i << "," << j<<","<<delt << endl;
                       // cout << "locali,localj:" << locali << "," << localj << endl;
                        if (locali < localj) {
                            no_tabu_pair.node1 = i;
                            no_tabu_pair.local1 = locali;
                            no_tabu_pair.node2 = j;
                            no_tabu_pair.local2 = localj;
                        } else {
                            no_tabu_pair.node1 = j;
                            no_tabu_pair.local1 = localj;
                            no_tabu_pair.node2 = i;
                            no_tabu_pair.local2 = locali;
                        }
                        no_tabu_samenumber = 1;
                        no_tabu_pair.delt = delt;
                    } else if (delt == no_tabu_pair.delt) {
                        tabu_samenumber++;
                        if (rand() % no_tabu_samenumber == 0) {
                            if (locali < localj) {
                                no_tabu_pair.node1 = i;
                                no_tabu_pair.local1 = locali;
                                no_tabu_pair.node2 = j;
                                no_tabu_pair.local2 = localj;
                            } else {
                                no_tabu_pair.node1 = j;
                                no_tabu_pair.local1 = localj;
                                no_tabu_pair.node2 = i;
                                no_tabu_pair.local2 = locali;
                            }
                        }
                    }
                }
            }
        }
    }
    if (tabu_pair.delt < no_tabu_pair.delt && tabu_pair.delt+local_solution < best_solution) {
        Pair = tabu_pair;
        cout << "** ";
    } else{
        Pair = no_tabu_pair;
    } //else return 0;//找不到可以改善最优解的交换节点对

    return 1;

}


void Solver::change_pair(Solution &sln, pair &Pair,ID nodeNum) {
    //更新tabu表
    TabuTenure[Pair.node2][Pair.node1] = TabuTenure[Pair.node1][Pair.node2] = 10 + iter + rand() % iter;
    vector <int> tempSolution_path = Solution_path;
    //cout << "node1:" << Pair.node1 << ",node2:" << Pair.node2 <<",delt:"<<Pair.delt<< endl;
        //顺序标记更新
    inpath[Solution_path[Pair.local1 - 1]][Pair.node1] = inpath[Pair.node1][Solution_path[Pair.local1 - 1]] = 0;
    inpath[Pair.node2][Solution_path[Pair.local2 + 1]] = inpath[Solution_path[Pair.local2 + 1]][Pair.node2] = 0;
    inpath[Solution_path[Pair.local1 - 1]][Pair.node2] = inpath[Pair.node2][Solution_path[Pair.local1 - 1]] = 1;
    inpath[Pair.node1][Solution_path[Pair.local2 + 1]] = inpath[Solution_path[Pair.local2 + 1]][Pair.node1] = 1;
    for (int i = Pair.local1, j = 0; i != Pair.local2 + 1; i++, j++) {
        Solution_path[i] = tempSolution_path[Pair.local2 - j];
    }
    local_solution += Pair.delt;
}

void Solver::check(Solution &sln,ID nodeNum) {
//   // cout << best_solution << endl;
    int best_cost = 0;
    for (int i = 1; i != nodeNum + 1; i++) {
        best_cost += aux.adjMat[Solution_path[i - 1]][Solution_path[i]];
    }
    if (best_cost != best_solution) {
        cout << "error" << endl; 
        cout << best_cost << endl;
    }
}
#pragma endregion Solver

}
