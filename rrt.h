#include <vector>
#include <queue>
#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#define PI 3.141592654

std::vector<double> randomNDimensionalVector(int n)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distribution(0.0, 2.0 * M_PI);

    std::vector<double> randomVector(n);

    for (int i = 0; i < n; ++i)
    {
        randomVector[i] = distribution(gen);
    }

    return randomVector;
}

float get_distance(std::vector<double> a, std::vector<double> b)
{
    float distance = 0.0;
    for (int i = 0; i < a.size(); ++i)
    {
        distance += pow(a[i] - b[i], 2);
    }
    return sqrt(distance);
}
struct cell
{
    std::vector<double> coordinates;
    cell *parent = nullptr;
};

class rrtplanner
{
public:
    std::vector<cell *> tree;
    int dimensions;
    int start_x, start_y, goal_x, goal_y;
    int x_size, y_size;
    int num_points{5000};
    double epsilon{0.5};
    double *map;
    bool gl{false};
    bool reached{false};
    std::vector<double> myloc;
    std::vector<double> goal;

    void generate_tree();
    void extend_tree(std::vector<double> randomVector, int loc);
    int find_nearest_neighbour(std::vector<double> randomVector);
    void backrtack(double ***plan, int *planlength);
    

    rrtplanner(int x_size, int y_size, int goal_x, int goal_y, int d, double *map, double *armstart_anglesV_rad,double *armgoal_anglesV_rad);
};


class rrtconnect
{
public:
    std::vector<cell *> start_tree,goal_tree;
    int dimensions;
    int start_x, start_y, goal_x, goal_y;
    int x_size, y_size;
    int num_points{10000};
    double epsilon{0.5};
    double *map;
    std::vector<double> myloc;
    std::vector<double> goal;
    bool start{true};
    bool trapped{false};
    bool reached{false};
    bool reverse{false};

    void generate_tree();
    void extend_tree_start(std::vector<double> randomVector, int loc);
    void extend_tree_end(std::vector<double> randomVector, int loc);
    void generatepath();
    int find_nearest_neighbour_start(std::vector<double> randomVector);
    int find_nearest_neighbour_end(std::vector<double> randomVector);
    void backrtack(double ***plan, int *planlength);
    void connect_start();
    void connect_end();
    

    rrtconnect(int x_size, int y_size, int goal_x, int goal_y, int d, double *map, double *armstart_anglesV_rad,double *armgoal_anglesV_rad);
};



class rrtstar
{
public:
    std::vector<cell *> tree;
    int dimensions;
    int start_x, start_y, goal_x, goal_y;
    int x_size, y_size;
    int num_points{5000};
    double epsilon{0.5};
    float prob{0.5};
    double *map;
    bool gl{false};
    bool reached{false};
    std::vector<double> myloc;
    std::vector<double> goal;

    void generate_tree();
    void extend_tree(std::vector<double> randomVector, int loc);
    int find_nearest_neighbour(std::vector<double> randomVector);
    void backrtack(double ***plan, int *planlength);
    

    rrtstar(int x_size, int y_size, int goal_x, int goal_y, int d, double *map, double *armstart_anglesV_rad,double *armgoal_anglesV_rad);
};