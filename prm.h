#include <vector>
#include <queue>
#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include <queue>
#define PI 3.141592654

struct pcell
{
    float h{1111110.0}, g{111111110.0}, f{111111110.0};
    // int parent_x, parent_y;
    std::vector<double> coordinates;
    // int c_x, c_y;
    bool collision{false}, visited{false};
    pcell *parent = nullptr;
    // int path_cost{0};
    std::vector<pcell *> neighbours;
};

class prmplanner
{
    public:
    std::vector<pcell *> tree;
    int numsamples{10000};

    int dimensions;
    int start_x, start_y, goal_x, goal_y;
    int x_size, y_size;
    std::vector<double> myloc;
    std::vector<double> goal;
    double *map;

    void generatetree();
    void connect(pcell * p);
    std::pair<int, int> neighbors();
    void catch_t(double ***plan, int *planlength);

    prmplanner(int x_size, int y_size, int goal_x, int goal_y, int d, double *map, double *armstart_anglesV_rad, double *armgoal_anglesV_rad);
};