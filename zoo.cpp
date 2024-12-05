//3E33912F8BAA7542FC4A1585D2DB6FE0312725B9

#include <iostream>
#include <vector>
#include <unordered_map>
#include <deque>
#include <getopt.h>
#include <limits>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace::std;

struct Cage {
    double x, y;
    int parent;    
    double distance = numeric_limits<double>::infinity();
    char type;
    bool visit = false;
};

double check_distance(const Cage &c1, const Cage &c2) {
    if (c1.type == c2.type || (c1.type == 'W' || c2.type == 'W')) {
        double e1 = c1.x - c2.x;
        double e2 = c1.y - c2.y;
        return e1 * e1 + e2 * e2;
    }
    return numeric_limits<double>::infinity();
}

double check_distance1(const pair<double, double> &c1, const pair<double, double> &c2) {
    double e1 = c1.first - c2.first;
    double e2 = c1.second - c2.second;
    return e1 * e1 + e2 * e2;
}

double check_distance2(const pair<double, double> &c1, const pair<double, double> &c2) {
    double e1 = c1.first - c2.first;
    double e2 = c1.second - c2.second;
    return sqrt(e1 * e1 + e2 * e2);
}

class ZooMSP {
 private:
    vector<Cage> CageList;
    double Weight = 0;

 public:
    void read_vertex();
    void make_mst();
    void print_mst();
};

void ZooMSP::read_vertex() {
    double xin, yin;
    size_t size;
    cin >> size;
    CageList.reserve(size);
    cin >> xin >> yin;
    while(!cin.fail()) {
        Cage c;
        c.x = xin;
        c.y = yin;
        if (xin > 0 || yin > 0) {
            c.type = 's';
        }
        else if(xin < 0 && yin < 0) {
            c.type = 'w';
        }
        else {
            c.type = 'W';
        }
        CageList.push_back(c);
        cin >> xin >> yin;
    }
}

void ZooMSP::make_mst() {
    CageList[0].distance = 0;
    int size = static_cast<int>(CageList.size());
    int count = 0;
    while (count < size) {
        double minDis = numeric_limits<double>::max();
        int min = -1;
        for (int i = 0; i < size; ++i) {
            if (!CageList[i].visit) {
                if (CageList[i].distance < minDis) {
                    minDis = CageList[i].distance;
                    min = static_cast<int>(i);
                }
            }
        }
        if (min == -1) {
            cerr << "Cannot construct MST\n";
            exit(1);
        }
        CageList[min].visit = true;
        for (int j = 0; j < size; ++j) {
            if (!CageList[j].visit) {
                double dis = check_distance(CageList[min], CageList[j]);
                if (dis < CageList[j].distance) {
                    CageList[j].distance = dis;
                    CageList[j].parent = min;
                }
            }
        }
        Weight += sqrt(CageList[min].distance);
        ++count;
    }
}

void ZooMSP::print_mst() {
    cout << Weight << "\n";
    int size = static_cast<int>(CageList.size());
    for (int i = 1; i < size; ++i) {
        if (i < CageList[i].parent) {
            cout << i << " " << CageList[i].parent << "\n";
        }
        else {
            cout << CageList[i].parent << " " << i << "\n";
        }
    }
}

class ZooTSP {
 private:
    int size;
    vector<pair<double, double>> CageList;
    vector<int> path;
    vector<int> currpath;
    vector<double> distance;
    vector<bool> visit;
    double curr_length = 0;
    double best_length = 0;
    bool promising(size_t permlength);
    double callmst(size_t permlength);

 public:
    void read_vertex();
    void make_tsp();
    void print_tsp();
    void setp();
    void genPerms(size_t permlength);
};

void ZooTSP::read_vertex() {
    double xin, yin;
    cin >> size;
    CageList.reserve(size);
    path.reserve(size + 1);
    cin >> xin >> yin;
    while(!cin.fail()) {
        CageList.push_back({xin, yin});
        cin >> xin >> yin;
    }
}

void ZooTSP::make_tsp() {
    path.push_back(0);
    path.push_back(1);
    path.push_back(0);
    size_t index_insert = 0;
    for (int j = 2; j < size; ++j) {
        double minCost = numeric_limits<double>::infinity();
        for (size_t i = 0; i < path.size() - 1; ++i) {
            double cost = check_distance2(CageList[path[i]], CageList[j]) + 
                          check_distance2(CageList[path[i + 1]], CageList[j]) -
                          check_distance2(CageList[path[i]], CageList[path[i + 1]]);
            if (cost < minCost) {
                minCost = cost;
                index_insert = i + 1;
            }
        }
        path.insert(path.begin() + index_insert, j);
    }
    path.pop_back();
}

void ZooTSP::print_tsp() {
    double total_distance = 0;
    for (int i = 0; i < size - 1; ++i) {
        total_distance += check_distance2(CageList[path[i]], CageList[path[i + 1]]);
    }
    total_distance += check_distance2(CageList[path[0]], CageList[path.back()]);
    cout << total_distance << "\n";
    for (int i = 0; i < size; ++i) {
        cout << path[i] << " ";
    }
    cout << "\n";
}

bool ZooTSP::promising(size_t permlength) {
    if (path.size() - permlength < 5)
        return true;

    if (curr_length > best_length)
        return false;

    double startbest = INFINITY;
    double endbest = INFINITY;
    size_t end = permlength - 1;
    for (size_t i = permlength; i < path.size(); ++i) {
        double cur = check_distance1(CageList[currpath[i]], CageList[currpath[0]]);
        if (cur < startbest)
            startbest = cur;

        cur = check_distance1(CageList[currpath[i]], CageList[currpath[end]]);;
        if (cur < endbest)
            endbest = cur;
    }
    startbest = sqrt(startbest);
    endbest = sqrt(endbest);
    double l = callmst(permlength);
    return curr_length + l + startbest + endbest < best_length;
}

void ZooTSP::genPerms(size_t permlength) {
    if (permlength == currpath.size()) {
        // Do something with the path
        double l = curr_length + check_distance2(CageList[path.front()], CageList[currpath.back()]);
        if (l < best_length) {
            best_length = l;
            path = currpath;
        }
        return;
    }
    if (!promising(permlength))
        return;
    for (size_t i = permlength; i < currpath.size(); ++i) {
        swap(currpath[permlength], currpath[i]);
        double l = check_distance2(CageList[currpath[permlength - 1]], CageList[currpath[permlength]]);
        curr_length += l;
        genPerms(permlength + 1);
        curr_length -= l;
        swap(currpath[permlength], currpath[i]);
    }
}

double ZooTSP::callmst(size_t permlength) {
    double weight = 0;
    size_t lefto = currpath.size() - permlength;
    distance.clear();
    visit.clear();
    distance.resize(lefto, numeric_limits<double>::infinity());
    visit.resize(lefto, false);
    distance[0] = 0;
    size_t min = 0;
    for (size_t count = 0; count < lefto; ++count) {
        double minDis = numeric_limits<double>::infinity();
        size_t nextmin;
        for (size_t i = 0; i < lefto; ++i) {
            if (!visit[i]) {
                double dis = check_distance1(CageList[currpath[min + permlength]], CageList[currpath[i + permlength]]);
                if (dis < distance[i])
                    distance[i] = dis;
                if (distance[i] < minDis) {
                    minDis = distance[i];
                    nextmin = i;
                }
            }
        }
        min = nextmin;
        visit[min] = true;
        weight += sqrt(distance[min]);
    }
    return weight;
}

void ZooTSP::setp() {
    currpath = path;
    for (int i = 0; i < size - 1; ++i) {
        best_length += check_distance2(CageList[path[i]], CageList[path[i + 1]]);
    }
    best_length += check_distance2(CageList[path[0]], CageList[path.back()]); 
}

int main(int argc, char **argv) {
    ios_base::sync_with_stdio(false);
    cout << setprecision(2);
    cout << fixed;
    option longOpts[] = {{"help", no_argument, nullptr, 'h'},
                        {"mode", required_argument, nullptr, 'm'},
                        { nullptr, 0, nullptr, '\0' }};
    int choice = 0, option_index = 0;
    string arg;
    while ((choice = getopt_long(argc, argv, "hm:", longOpts, &option_index)) != -1) {
        switch (choice) {
            case 'm':
                arg = optarg;
                if (arg == "MST") {
                    ZooMSP zoo;
                    zoo.read_vertex();
                    zoo.make_mst();
                    zoo.print_mst();
                }
                else if (arg == "FASTTSP") {
                    ZooTSP zoo;
                    zoo.read_vertex();
                    zoo.make_tsp();
                    zoo.print_tsp();
                }
                else if (arg == "OPTTSP") {
                    ZooTSP zoo;
                    zoo.read_vertex();
                    zoo.make_tsp();
                    zoo.setp();
                    zoo.genPerms(1);
                    zoo.print_tsp();
                }
                else {
                    cerr << "Invalid argument for mode.\n";
                    exit(1);
                }
                break;
            case 'h':
                cout << "This is zoo.\n";
                exit(0);
            default:
                cerr << "invalid command line argument\n";
                exit(1);
        }
    }
}
