/*=================================================================
 *
 * planner.c
 *
 *=================================================================*/
#include <math.h>
#include <random>
#include <vector>
#include <array>
#include <algorithm>
#include <chrono>
#include <thread>

#include <tuple>
#include <string>
#include <stdexcept>
#include <regex>	// For regex and split logic
#include <iostream> // cout, endl
#include <fstream>	// For reading/writing files
#include <assert.h>
#include "rrt.h"
#include "prm.h"
/* Input Arguments */
#define MAP_IN prhs[0]
#define ARMSTART_IN prhs[1]
#define ARMGOAL_IN prhs[2]
#define PLANNER_ID_IN prhs[3]

/* Planner Ids */
#define RRT 0
#define RRTCONNECT 1
#define RRTSTAR 2
#define PRM 3

/* Output Arguments */
#define PLAN_OUT plhs[0]
#define PLANLENGTH_OUT plhs[1]

#define GETMAPINDEX(X, Y, XSIZE, YSIZE) (Y * XSIZE + X)

#if !defined(MAX)
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif

#define PI 3.141592654

// the length of each link in the arm
#define LINKLENGTH_CELLS 10

// Some potentially helpful imports
using std::array;
using std::cout;
using std::endl;
using std::make_tuple;
using std::runtime_error;
using std::string;
using std::tie;
using std::tuple;
using std::vector;

struct CompareCellByF
{
	bool operator()(const pcell *a, const pcell *b)
	{
		return a->f > b->f; // Change to '<' for min priority queue
	}
};

//*******************************************************************************************************************//
//                                                                                                                   //
//                                                GIVEN FUNCTIONS                                                    //
//                                                                                                                   //
//*******************************************************************************************************************//

/// @brief
/// @param filepath
/// @return map, x_size, y_size
tuple<double *, int, int> loadMap(string filepath)
{
	std::FILE *f = fopen(filepath.c_str(), "r");
	if (f)
	{
	}
	else
	{
		printf("Opening file failed! \n");
		throw runtime_error("Opening map file failed!");
	}
	int height, width;
	if (fscanf(f, "height %d\nwidth %d\n", &height, &width) != 2)
	{
		throw runtime_error("Invalid loadMap parsing map metadata");
	}

	////// Go through file and add to m_occupancy
	double *map = new double[height * width];

	double cx, cy, cz;
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			char c;
			do
			{
				if (fscanf(f, "%c", &c) != 1)
				{
					throw runtime_error("Invalid parsing individual map data");
				}
			} while (isspace(c));
			if (!(c == '0'))
			{
				map[y + x * width] = 1; // Note transposed from visual
			}
			else
			{
				map[y + x * width] = 0;
			}
		}
	}
	fclose(f);
	return make_tuple(map, width, height);
}

// Splits string based on deliminator
vector<string> split(const string &str, const string &delim)
{
	// https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c/64886763#64886763
	const std::regex ws_re(delim);
	return {std::sregex_token_iterator(str.begin(), str.end(), ws_re, -1), std::sregex_token_iterator()};
}

double *doubleArrayFromString(string str)
{
	vector<string> vals = split(str, ",");
	double *ans = new double[vals.size()];
	for (int i = 0; i < vals.size(); ++i)
	{
		ans[i] = std::stod(vals[i]);
	}
	return ans;
}

bool equalDoubleArrays(double *v1, double *v2, int size)
{
	for (int i = 0; i < size; ++i)
	{
		if (abs(v1[i] - v2[i]) > 1e-3)
		{
			cout << endl;
			return false;
		}
	}
	return true;
}

typedef struct
{
	int X1, Y1;
	int X2, Y2;
	int Increment;
	int UsingYIndex;
	int DeltaX, DeltaY;
	int DTerm;
	int IncrE, IncrNE;
	int XIndex, YIndex;
	int Flipped;
} bresenham_param_t;

void ContXY2Cell(double x, double y, short unsigned int *pX, short unsigned int *pY, int x_size, int y_size)
{
	double cellsize = 1.0;
	// take the nearest cell
	*pX = (int)(x / (double)(cellsize));
	if (x < 0)
		*pX = 0;
	if (*pX >= x_size)
		*pX = x_size - 1;

	*pY = (int)(y / (double)(cellsize));
	if (y < 0)
		*pY = 0;
	if (*pY >= y_size)
		*pY = y_size - 1;
}

void get_bresenham_parameters(int p1x, int p1y, int p2x, int p2y, bresenham_param_t *params)
{
	params->UsingYIndex = 0;

	if (fabs((double)(p2y - p1y) / (double)(p2x - p1x)) > 1)
		(params->UsingYIndex)++;

	if (params->UsingYIndex)
	{
		params->Y1 = p1x;
		params->X1 = p1y;
		params->Y2 = p2x;
		params->X2 = p2y;
	}
	else
	{
		params->X1 = p1x;
		params->Y1 = p1y;
		params->X2 = p2x;
		params->Y2 = p2y;
	}

	if ((p2x - p1x) * (p2y - p1y) < 0)
	{
		params->Flipped = 1;
		params->Y1 = -params->Y1;
		params->Y2 = -params->Y2;
	}
	else
		params->Flipped = 0;

	if (params->X2 > params->X1)
		params->Increment = 1;
	else
		params->Increment = -1;

	params->DeltaX = params->X2 - params->X1;
	params->DeltaY = params->Y2 - params->Y1;

	params->IncrE = 2 * params->DeltaY * params->Increment;
	params->IncrNE = 2 * (params->DeltaY - params->DeltaX) * params->Increment;
	params->DTerm = (2 * params->DeltaY - params->DeltaX) * params->Increment;

	params->XIndex = params->X1;
	params->YIndex = params->Y1;
}

void get_current_point(bresenham_param_t *params, int *x, int *y)
{
	if (params->UsingYIndex)
	{
		*y = params->XIndex;
		*x = params->YIndex;
		if (params->Flipped)
			*x = -*x;
	}
	else
	{
		*x = params->XIndex;
		*y = params->YIndex;
		if (params->Flipped)
			*y = -*y;
	}
}

int get_next_point(bresenham_param_t *params)
{
	if (params->XIndex == params->X2)
	{
		return 0;
	}
	params->XIndex += params->Increment;
	if (params->DTerm < 0 || (params->Increment < 0 && params->DTerm <= 0))
		params->DTerm += params->IncrE;
	else
	{
		params->DTerm += params->IncrNE;
		params->YIndex += params->Increment;
	}
	return 1;
}

int IsValidLineSegment(double x0, double y0, double x1, double y1, double *map,
					   int x_size, int y_size)
{
	bresenham_param_t params;
	int nX, nY;
	short unsigned int nX0, nY0, nX1, nY1;

	// printf("checking link <%f %f> to <%f %f>\n", x0,y0,x1,y1);

	// make sure the line segment is inside the environment
	if (x0 < 0 || x0 >= x_size ||
		x1 < 0 || x1 >= x_size ||
		y0 < 0 || y0 >= y_size ||
		y1 < 0 || y1 >= y_size)
		return 0;

	ContXY2Cell(x0, y0, &nX0, &nY0, x_size, y_size);
	ContXY2Cell(x1, y1, &nX1, &nY1, x_size, y_size);

	// printf("checking link <%d %d> to <%d %d>\n", nX0,nY0,nX1,nY1);

	// iterate through the points on the segment
	get_bresenham_parameters(nX0, nY0, nX1, nY1, &params);
	do
	{
		get_current_point(&params, &nX, &nY);
		if (map[GETMAPINDEX(nX, nY, x_size, y_size)] == 1)
			return 0;
	} while (get_next_point(&params));

	return 1;
}

int IsValidArmConfiguration(double *angles, int numofDOFs, double *map,
							int x_size, int y_size)
{
	double x0, y0, x1, y1;
	int i;

	// iterate through all the links starting with the base
	x1 = ((double)x_size) / 2.0;
	y1 = 0;
	for (i = 0; i < numofDOFs; i++)
	{
		// compute the corresponding line segment
		x0 = x1;
		y0 = y1;
		x1 = x0 + LINKLENGTH_CELLS * cos(2 * PI - angles[i]);
		y1 = y0 - LINKLENGTH_CELLS * sin(2 * PI - angles[i]);

		// check the validity of the corresponding line segment
		if (!IsValidLineSegment(x0, y0, x1, y1, map, x_size, y_size))
			return 0;
	}
	return 1;
}

//*******************************************************************************************************************//
//                                                                                                                   //
//                                          DEFAULT PLANNER FUNCTION                                                 //
//                                                                                                                   //
//*******************************************************************************************************************//

static void planner(
	double *map,
	int x_size,
	int y_size,
	double *armstart_anglesV_rad,
	double *armgoal_anglesV_rad,
	int numofDOFs,
	double ***plan,
	int *planlength)
{
	// no plan by default
	*plan = NULL;
	*planlength = 0;

	// for now just do straight interpolation between start and goal checking for the validity of samples

	double distance = 0;
	int i, j;
	for (j = 0; j < numofDOFs; j++)
	{
		if (distance < fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]))
			distance = fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]);
	}
	int numofsamples = (int)(distance / (PI / 20));
	if (numofsamples < 2)
	{
		printf("the arm is already at the goal\n");
		return;
	}
	*plan = (double **)malloc(numofsamples * sizeof(double *));
	int firstinvalidconf = 1;
	for (i = 0; i < numofsamples; i++)
	{
		(*plan)[i] = (double *)malloc(numofDOFs * sizeof(double));
		for (j = 0; j < numofDOFs; j++)
		{
			(*plan)[i][j] = armstart_anglesV_rad[j] + ((double)(i) / (numofsamples - 1)) * (armgoal_anglesV_rad[j] - armstart_anglesV_rad[j]);
		}
		if (!IsValidArmConfiguration((*plan)[i], numofDOFs, map, x_size, y_size) && firstinvalidconf)
		{
			firstinvalidconf = 1;
			printf("ERROR: Invalid arm configuration!!!\n");
		}
	}
	*planlength = numofsamples;

	return;
}

//*******************************************************************************************************************//
//                                                                                                                   //
//                                              RRT IMPLEMENTATION                                                   //
//                                                                                                                   //
//*******************************************************************************************************************//

rrtplanner::rrtplanner(int x_size, int y_size, int goal_x, int goal_y, int d, double *map, double *armstart_anglesV_rad, double *armgoal_anglesV_rad)
{
	this->x_size = x_size;
	this->y_size = y_size;
	this->goal_x = goal_x;
	this->goal_y = goal_y;
	this->dimensions = d;
	this->map = map;
	for (int i = 0; i < dimensions; i++)
	{
		this->myloc.push_back(armstart_anglesV_rad[i]);
		this->goal.push_back(armgoal_anglesV_rad[i]);
	}
}

void rrtplanner::generate_tree()
{

	std::vector<double> randomVector(dimensions);
	int loc{0};
	cell *start = new cell;
	start->coordinates = myloc;
	tree.push_back(start);
	for (int i = 0; i < num_points; i++)
	{
		randomVector = randomNDimensionalVector(dimensions);
		loc = find_nearest_neighbour(randomVector);
		extend_tree(randomVector, loc);
		if (reached)
			break;
	}
}
int rrtplanner::find_nearest_neighbour(std::vector<double> randomVector)
{
	float min_dist{100000000.0};
	int loc{0};
	float dist{0};
	for (int i = 0; i < tree.size(); i++)
	{
		dist = get_distance(tree[i]->coordinates, randomVector);
		if (dist < min_dist)
		{
			min_dist = dist;
			loc = i;
		}
	}
	return loc;
}

void rrtplanner::extend_tree(std::vector<double> randomVector, int loc)
{
	cell *c = new cell;
	c->coordinates = randomVector;
	std::vector<double> interpol(dimensions);
	std::vector<double> dirvec(dimensions);
	std::vector<double> safe(dimensions);
	safe = tree[loc]->coordinates;
	bool flag{true};
	int numofsamples;
	double distance = 0;
	double *pointer;
	int i, j;
	for (j = 0; j < dimensions; j++)
	{
		if (distance < fabs(tree[loc]->coordinates[j] - randomVector[j]))
		{
			distance = fabs(tree[loc]->coordinates[j] - randomVector[j]);
		}
		if ((randomVector[j] - tree[loc]->coordinates[j]) < 0)
		{
			dirvec[j] = std::max((randomVector[j] - tree[loc]->coordinates[j]), -epsilon);
		}
		else
		{
			dirvec[j] = std::min((randomVector[j] - tree[loc]->coordinates[j]), epsilon);
		}
	}
	numofsamples = (int)(distance / (PI / 10));
	if (numofsamples < 2)
	{
		// printf("the arm is already at the goal\n");
		return;
		flag = false;
	}

	for (i = 0; i < numofsamples; i++)
	{

		for (j = 0; j < dimensions; j++)
		{
			interpol[j] = tree[loc]->coordinates[j] + ((double)(i) / (numofsamples - 1)) * dirvec[j];
		}
		pointer = interpol.data();
		if (!IsValidArmConfiguration(pointer, dimensions, map, x_size, y_size))
		{
			break;
			printf("ERROR: Invalid arm configuration!!!\n");
		}
		safe = interpol;
	}
	c->coordinates = safe;
	c->parent = tree[loc];
	tree.push_back(c);
	if (get_distance(c->coordinates, goal) < epsilon * 2)
	{
		reached = true;
	}
}

void rrtplanner::backrtack(double ***plan, int *planlength)
{
	std::vector<double> goal(dimensions);
	goal = this->goal;
	int loc = find_nearest_neighbour(goal);
	float dist = get_distance(tree[loc]->coordinates, goal);
	if (dist > 2.0)
	{
		std::cout << "path not found" << std::endl;
		return;
	}
	cell *c = tree[loc];
	cell *g = new cell;
	g->coordinates = goal;
	g->parent = c;

	std::vector<cell *> path;
	int i{0};
	while (g->parent != nullptr)
	{
		path.push_back(g);
		g = g->parent;
		i++;
	}
	path.push_back(g);
	std::reverse(path.begin(), path.end());

	*plan = (double **)malloc(path.size() * sizeof(double *));
	int j;
	for (j = 0; j < path.size(); j++)
	{
		(*plan)[j] = (double *)malloc(dimensions * sizeof(double));
		(*plan)[j] = path[j]->coordinates.data();
	}
	*planlength = path.size();
}

static void plannerRRT(
	double *map,
	int x_size,
	int y_size,
	double *armstart_anglesV_rad,
	double *armgoal_anglesV_rad,
	int numofDOFs,
	double ***plan,
	int *planlength)
{
	/* TODO: Replace with your implementation */
	// planner(map, x_size, y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, plan, planlength);
	rrtplanner rrt(x_size, y_size, armgoal_anglesV_rad[0], armgoal_anglesV_rad[1], numofDOFs, map, armstart_anglesV_rad, armgoal_anglesV_rad);
	rrt.generate_tree();
	rrt.backrtack(plan, planlength);
}

//*******************************************************************************************************************//
//                                                                                                                   //
//                                         RRT CONNECT IMPLEMENTATION                                                //
//                                                                                                                   //
//*******************************************************************************************************************//
rrtconnect::rrtconnect(int x_size, int y_size, int goal_x, int goal_y, int d, double *map, double *armstart_anglesV_rad, double *armgoal_anglesV_rad)
{
	this->x_size = x_size;
	this->y_size = y_size;
	this->goal_x = goal_x;
	this->goal_y = goal_y;
	this->dimensions = d;
	this->map = map;
	for (int i = 0; i < dimensions; i++)
	{
		this->myloc.push_back(armstart_anglesV_rad[i]);
		this->goal.push_back(armgoal_anglesV_rad[i]);
	}
}

void rrtconnect::generate_tree()
{

	std::vector<double> randomVector(dimensions);
	int loc{0};
	cell *st = new cell;
	st->coordinates = myloc;
	start_tree.push_back(st);
	cell *end = new cell;
	end->coordinates = goal;
	goal_tree.push_back(end);
	for (int i = 0; i < num_points; i++)
	{
		randomVector = randomNDimensionalVector(dimensions);

		loc = find_nearest_neighbour_start(randomVector);
		extend_tree_start(randomVector, loc);
		connect_start();

		if (reached)
			break;

		loc = find_nearest_neighbour_end(randomVector);
		extend_tree_end(randomVector, loc);
		connect_end();
		if (reached)
			break;
	}
}
void rrtconnect::connect_start()
{
	int loc{0};
	int prev_loc{-1};
	std::vector<double> coor(dimensions);
	coor = start_tree.back()->coordinates;
	loc = find_nearest_neighbour_end(coor);
	while (true)
	{
		trapped = false;
		extend_tree_end(coor, loc);
		loc = find_nearest_neighbour_end(coor);
		if ((prev_loc == loc) || trapped)
		{
			if (prev_loc == loc)
				reached = true;
			break;
		}
		prev_loc = loc;
	}
}
void rrtconnect::connect_end()
{
	int loc{0};
	int prev_loc{-1};
	std::vector<double> coor(dimensions);
	coor = goal_tree.back()->coordinates;
	loc = find_nearest_neighbour_start(coor);
	while (true)
	{
		trapped = false;
		extend_tree_start(coor, loc);
		loc = find_nearest_neighbour_start(coor);
		if ((prev_loc == loc) || trapped)
		{
			if (prev_loc == loc)
				reached = true;
			break;
		}
		prev_loc = loc;
	}
}
int rrtconnect::find_nearest_neighbour_start(std::vector<double> randomVector)
{
	float min_dist{100000000.0};
	int loc{0};
	for (int i = 0; i < start_tree.size(); i++)
	{
		float dist = get_distance(start_tree[i]->coordinates, randomVector);
		if (dist < min_dist)
		{
			min_dist = dist;
			loc = i;
		}
	}
	return loc;
}
int rrtconnect::find_nearest_neighbour_end(std::vector<double> randomVector)
{
	float min_dist{100000000.0};
	int loc{0};
	for (int i = 0; i < goal_tree.size(); i++)
	{
		float dist = get_distance(goal_tree[i]->coordinates, randomVector);
		if (dist < min_dist)
		{
			min_dist = dist;
			loc = i;
		}
	}
	return loc;
}

void rrtconnect::extend_tree_start(std::vector<double> randomVector, int loc)
{
	cell *c = new cell;
	c->coordinates = randomVector;
	std::vector<double> interpol(dimensions);
	std::vector<double> dirvec(dimensions);
	std::vector<double> safe(dimensions);
	safe = start_tree[loc]->coordinates;
	bool flag{true};
	double distance = 0;
	int i, j;
	for (j = 0; j < dimensions; j++)
	{
		if (distance < fabs(start_tree[loc]->coordinates[j] - randomVector[j]))
		{
			distance = fabs(start_tree[loc]->coordinates[j] - randomVector[j]);
		}
		if ((randomVector[j] - start_tree[loc]->coordinates[j]) < 0)
		{
			dirvec[j] = std::max((randomVector[j] - start_tree[loc]->coordinates[j]), -epsilon);
		}
		else
		{
			dirvec[j] = std::min((randomVector[j] - start_tree[loc]->coordinates[j]), epsilon);
		}
	}
	int numofsamples = (int)(distance / (PI / 15));
	if (numofsamples < 2)
	{
		return;
	}

	for (i = 0; i < numofsamples; i++)
	{

		for (j = 0; j < dimensions; j++)
		{
			interpol[j] = start_tree[loc]->coordinates[j] + ((double)(i) / (numofsamples - 1)) * dirvec[j];
		}
		double *pointer = interpol.data();
		if (!IsValidArmConfiguration(pointer, dimensions, map, x_size, y_size))
		{
			trapped = true;
			break;
			printf("ERROR: Invalid arm configuration!!!\n");
		}
		safe = interpol;
	}
	c->coordinates = safe;
	c->parent = start_tree[loc];
	start_tree.push_back(c);
}

void rrtconnect::extend_tree_end(std::vector<double> randomVector, int loc)
{
	cell *c = new cell;
	c->coordinates = randomVector;
	std::vector<double> interpol(dimensions);
	std::vector<double> dirvec(dimensions);
	std::vector<double> safe(dimensions);
	safe = goal_tree[loc]->coordinates;
	bool flag{true};
	double distance = 0;
	int i, j;
	for (j = 0; j < dimensions; j++)
	{
		if (distance < fabs(goal_tree[loc]->coordinates[j] - randomVector[j]))
		{
			distance = fabs(goal_tree[loc]->coordinates[j] - randomVector[j]);
		}
		if ((randomVector[j] - goal_tree[loc]->coordinates[j]) < 0)
		{
			dirvec[j] = std::max((randomVector[j] - goal_tree[loc]->coordinates[j]), -epsilon);
		}
		else
		{
			dirvec[j] = std::min((randomVector[j] - goal_tree[loc]->coordinates[j]), epsilon);
		}
	}
	int numofsamples = (int)(distance / (PI / 15));
	if (numofsamples < 2)
	{
		// printf("the arm is already at the goal\n");
		return;
	}

	for (i = 0; i < numofsamples; i++)
	{

		for (j = 0; j < dimensions; j++)
		{
			interpol[j] = goal_tree[loc]->coordinates[j] + ((double)(i) / (numofsamples - 1)) * dirvec[j];
		}
		double *pointer = interpol.data();
		if (!IsValidArmConfiguration(pointer, dimensions, map, x_size, y_size))
		{
			trapped = true;

			break;
			printf("ERROR: Invalid arm configuration!!!\n");
		}
		safe = interpol;
	}
	c->coordinates = safe;
	c->parent = goal_tree[loc];
	goal_tree.push_back(c);
}

void rrtconnect::backrtack(double ***plan, int *planlength)
{
	std::vector<cell *> path;
	int i;
	cell *c = start_tree.back();
	cell *g = goal_tree.back();

	while (g->parent != nullptr)
	{
		path.push_back(g);
		g = g->parent;
		i++;
	}
	path.push_back(g);
	std::reverse(path.begin(), path.end());
	while (c->parent != nullptr)
	{
		path.push_back(c);
		c = c->parent;
		i++;
	}
	path.push_back(c);
	std::reverse(path.begin(), path.end());
	*plan = (double **)malloc(path.size() * sizeof(double *));
	int j;
	for (j = 0; j < path.size(); j++)
	{
		(*plan)[j] = (double *)malloc(dimensions * sizeof(double));
		(*plan)[j] = path[j]->coordinates.data();
	}
	*planlength = path.size();
}

static void plannerRRTConnect(
	double *map,
	int x_size,
	int y_size,
	double *armstart_anglesV_rad,
	double *armgoal_anglesV_rad,
	int numofDOFs,
	double ***plan,
	int *planlength)
{
	/* TODO: Replace with your implementation */
	// planner(map, x_size, y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, plan, planlength);
	rrtconnect rrc(x_size, y_size, armgoal_anglesV_rad[0], armgoal_anglesV_rad[1], numofDOFs, map, armstart_anglesV_rad, armgoal_anglesV_rad);
	rrc.generate_tree();
	rrc.backrtack(plan, planlength);
}

//*******************************************************************************************************************//
//                                                                                                                   //
//                                           RRT STAR IMPLEMENTATION                                                 //
//                                                                                                                   //
//*******************************************************************************************************************//
rrtstar::rrtstar(int x_size, int y_size, int goal_x, int goal_y, int d, double *map, double *armstart_anglesV_rad, double *armgoal_anglesV_rad)
{
	this->x_size = x_size;
	this->y_size = y_size;
	this->goal_x = goal_x;
	this->goal_y = goal_y;
	this->dimensions = d;
	this->map = map;
	for (int i = 0; i < dimensions; i++)
	{
		this->myloc.push_back(armstart_anglesV_rad[i]);
		this->goal.push_back(armgoal_anglesV_rad[i]);
	}
}

void rrtstar::generate_tree()
{

	std::vector<double> randomVector(dimensions);
	int loc{0};
	cell *start = new cell;
	start->coordinates = myloc;
	tree.push_back(start);
	for (int i = 0; i < num_points; i++)
	{
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<double> distribution(0.0, 1.0);
		if (distribution(gen) < prob)
		{
			randomVector = goal;
		}

		else
		{
			randomVector = randomNDimensionalVector(dimensions);
		}
		if (reached)
			break;
		loc = find_nearest_neighbour(randomVector);
		extend_tree(randomVector, loc);
	}
}
int rrtstar::find_nearest_neighbour(std::vector<double> randomVector)
{
	float min_dist{100000000.0};
	int loc{0};
	for (int i = 0; i < tree.size(); i++)
	{
		float dist = get_distance(tree[i]->coordinates, randomVector);
		if (dist < min_dist)
		{
			min_dist = dist;
			loc = i;
		}
	}
	return loc;
}

void rrtstar::extend_tree(std::vector<double> randomVector, int loc)
{
	cell *c = new cell;
	c->coordinates = randomVector;
	std::vector<double> interpol(dimensions);
	std::vector<double> dirvec(dimensions);
	std::vector<double> safe(dimensions);
	;
	safe = tree[loc]->coordinates;
	bool flag{false};
	double distance = 0;
	float dist ;
	int i, j;
	for (j = 0; j < dimensions; j++)
	{
		if (distance < fabs(tree[loc]->coordinates[j] - randomVector[j]))
		{
			distance = fabs(tree[loc]->coordinates[j] - randomVector[j]);
		}
		if ((randomVector[j] - tree[loc]->coordinates[j]) < 0)
		{
			dirvec[j] = std::max((randomVector[j] - tree[loc]->coordinates[j]), -epsilon);
		}
		else
		{
			dirvec[j] = std::min((randomVector[j] - tree[loc]->coordinates[j]), epsilon);
		}
	}
	int numofsamples = (int)(distance / (PI / 10));
	if (numofsamples < 2)
	{
		// printf("the arm is already at the goal\n");
		return;
	}

	for (i = 0; i < numofsamples; i++)
	{

		for (j = 0; j < dimensions; j++)
		{
			interpol[j] = tree[loc]->coordinates[j] + ((double)(i) / (numofsamples - 1)) * dirvec[j];
		}
		double *pointer = interpol.data();
		if (!IsValidArmConfiguration(pointer, dimensions, map, x_size, y_size))
		{
			flag = true;
			break;
			printf("ERROR: Invalid arm configuration!!!\n");
		}
		safe = interpol;
	}
	dist = get_distance(safe, goal);
	if (dist < 2*epsilon)
	{
		safe = goal;
		reached = true;
	}
	c->coordinates = safe;
	c->parent = tree[loc];
	tree.push_back(c);
}

void rrtstar::backrtack(double ***plan, int *planlength)
{
	std::vector<double> goal(dimensions);
	goal = this->goal;
	int loc = find_nearest_neighbour(goal);
	float dist = get_distance(tree[loc]->coordinates, goal);
	// cout << dist << "  dist" << endl;
	if (dist > 2.0)
	{
		std::cout << "path not found" << std::endl;
		return;
	}
	cell *c = tree[loc];
	cell *g = new cell;
	g->coordinates = goal;
	g->parent = c;

	std::vector<cell *> path;
	int i{0};
	while (g->parent != nullptr)
	{
		path.push_back(g);
		g = g->parent;
		i++;
	}
	path.push_back(g);
	std::reverse(path.begin(), path.end());

	*plan = (double **)malloc(path.size() * sizeof(double *));
	int j;
	for (j = 0; j < path.size(); j++)
	{
		(*plan)[j] = (double *)malloc(dimensions * sizeof(double));
		(*plan)[j] = path[j]->coordinates.data();
	}
	*planlength = path.size();
}

static void plannerRRTStar(
	double *map,
	int x_size,
	int y_size,
	double *armstart_anglesV_rad,
	double *armgoal_anglesV_rad,
	int numofDOFs,
	double ***plan,
	int *planlength)
{
	/* TODO: Replace with your implementation */
	// planner(map, x_size, y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, plan, planlength);
	rrtstar rrs(x_size, y_size, armgoal_anglesV_rad[0], armgoal_anglesV_rad[1], numofDOFs, map, armstart_anglesV_rad, armgoal_anglesV_rad);
	rrs.generate_tree();
	rrs.backrtack(plan, planlength);
}

//*******************************************************************************************************************//
//                                                                                                                   //
//                                              PRM IMPLEMENTATION                                                   //
//                                                                                                                   //
//*******************************************************************************************************************//

prmplanner::prmplanner(int x_size, int y_size, int goal_x, int goal_y, int d, double *map, double *armstart_anglesV_rad, double *armgoal_anglesV_rad)
{
	this->x_size = x_size;
	this->y_size = y_size;
	this->goal_x = goal_x;
	this->goal_y = goal_y;
	this->dimensions = d;
	this->map = map;
	for (int i = 0; i < dimensions; i++)
	{
		this->myloc.push_back(armstart_anglesV_rad[i]);
		this->goal.push_back(armgoal_anglesV_rad[i]);
	}
}

void prmplanner::generatetree()
{
	vector<double> randomVector(dimensions);
	double *pointer;
	for (int i = 0; i < numsamples; i++)
	{
		randomVector = randomNDimensionalVector(dimensions);

		pointer = randomVector.data();
		if (IsValidArmConfiguration(pointer, dimensions, map, x_size, y_size))
		{
			pcell *p = new pcell;
			p->coordinates = randomVector;
			connect(p);
		}
	}
}

void prmplanner::connect(pcell *p)
{
	float dist;
	std::vector<double> interpol(dimensions);
	bool flag{true};
	int i, j;
	double *pointer;
	double distance = 0;
	for (int i = 0; i < tree.size(); i++)
	{
		dist = get_distance(tree[i]->coordinates, p->coordinates);
		if (dist < 2)
		{
			for (j = 0; j < dimensions; j++)
			{
				if (distance < fabs(tree[i]->coordinates[j] - p->coordinates[j]))
				{
					distance = fabs(tree[i]->coordinates[j] - p->coordinates[j]);
				}
			}
			int numofsamples = (int)(distance / (PI / 20));
			if (numofsamples < 2)
			{
				// printf("the arm is already at the goal\n");
				return;
			}
			for (int k = 0; k < numofsamples; k++)
			{
				for (j = 0; j < dimensions; j++)
				{
					interpol[j] = tree[i]->coordinates[j] + ((double)(k) / (numofsamples - 1)) * (p->coordinates[j] - tree[i]->coordinates[j]);
				}
				pointer = interpol.data();
				if (!IsValidArmConfiguration(pointer, dimensions, map, x_size, y_size))
				{
					flag = false;
					break;
					printf("ERROR: Invalid arm configuration!!!\n");
				}
			}
			if (flag)
			{
				tree[i]->neighbours.push_back(p);
				p->neighbours.push_back(tree[i]);
			}
		}
		flag = true;
	}
	dist = get_distance(goal, p->coordinates);
	p->h = dist;

	tree.push_back(p);
}

void prmplanner::catch_t(double ***plan, int *planlength)
{
	pcell *s = new pcell;
	pcell *e = new pcell;
	double sd{100000}, se{1000000};
	s->coordinates = myloc;
	e->coordinates = goal;
	float dist;
	int sloc, eloc;
	for (int i = 0; i < tree.size(); i++)
	{
		dist = get_distance(tree[i]->coordinates, myloc);
		if (dist < sd)
		{
			sloc = i;
			sd = dist;
		}
		dist = get_distance(tree[i]->coordinates, goal);
		if (dist < se)
		{
			eloc = i;
			se = dist;
		}
	}
	s->neighbours.push_back(tree[sloc]);
	float gar = get_distance(tree[sloc]->coordinates, s->coordinates);
	gar = get_distance(tree[eloc]->coordinates, e->coordinates);
	tree[sloc]->neighbours.push_back(s);
	e->neighbours.push_back(tree[eloc]);
	tree[eloc]->neighbours.push_back(e);
	tree.push_back(s);
	tree.push_back(e);
	std::this_thread::sleep_for(std::chrono::milliseconds(1000));

	std::vector<double> var(dimensions);

	std::priority_queue<pcell *, vector<pcell *>, CompareCellByF> priorityQueue;
	s->g = 0;
	s->f = s->g + s->h;
	priorityQueue.push(s);

	pcell *current = s;
	vector<pcell *> path;

	while (priorityQueue.size() > 0)
	{
		current = priorityQueue.top();
		current->visited = true;
		priorityQueue.pop();
		if (current->coordinates == goal)
		{
			while (current->parent != nullptr)
			{
				path.push_back(current);
				current = current->parent;
			}
			path.push_back(current);
			break;
		}
		var = current->coordinates;

		for (int k = 0; k < current->neighbours.size(); k++)
		{
			dist = get_distance(current->coordinates, current->neighbours[k]->coordinates);
			if (!current->neighbours[k]->visited && current->neighbours[k]->g > (current->g + dist))
			{
				current->neighbours[k]->visited = true;
				current->neighbours[k]->g = current->g + dist;
				current->neighbours[k]->f = current->g + current->h;
				current->neighbours[k]->parent = current;
				priorityQueue.push(current->neighbours[k]);
			}
		}
	}
	std::reverse(path.begin(), path.end());
	*plan = (double **)malloc(path.size() * sizeof(double *));
	int j;
	for (j = 0; j < path.size(); j++)
	{
		(*plan)[j] = (double *)malloc(dimensions * sizeof(double));
		(*plan)[j] = path[j]->coordinates.data();
	}
	*planlength = path.size();
}

static void plannerPRM(
	double *map,
	int x_size,
	int y_size,
	double *armstart_anglesV_rad,
	double *armgoal_anglesV_rad,
	int numofDOFs,
	double ***plan,
	int *planlength)
{
	/* TODO: Replace with your implementation */
	prmplanner prm(x_size, y_size, armgoal_anglesV_rad[0], armgoal_anglesV_rad[1], numofDOFs, map, armstart_anglesV_rad, armgoal_anglesV_rad);
	prm.generatetree();
	prm.catch_t(plan, planlength);
	// planner(map, x_size, y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, plan, planlength);
}

//*******************************************************************************************************************//
//                                                                                                                   //
//                                                MAIN FUNCTION                                                      //
//                                                                                                                   //
//*******************************************************************************************************************//

/** Your final solution will be graded by an grading script which will
 * send the default 6 arguments:
 *    map, numOfDOFs, commaSeparatedStartPos, commaSeparatedGoalPos,
 *    whichPlanner, outputFilePath
 * An example run after compiling and getting the planner.out executable
 * >> ./planner.out map1.txt 5 1.57,0.78,1.57,0.78,1.57 0.392,2.35,3.14,2.82,4.71 0 output.txt
 * See the hw handout for full information.
 * If you modify this for testing (e.g. to try out different hyper-parameters),
 * make sure it can run with the original 6 commands.
 * Programs that do not will automatically get a 0.
 * */
int main(int argc, char **argv)
{
	double *map;
	int x_size, y_size;

	tie(map, x_size, y_size) = loadMap(argv[1]);
	const int numOfDOFs = std::stoi(argv[2]);
	double *startPos = doubleArrayFromString(argv[3]);
	double *goalPos = doubleArrayFromString(argv[4]);
	int whichPlanner = std::stoi(argv[5]);
	string outputFile = argv[6];

	if (!IsValidArmConfiguration(startPos, numOfDOFs, map, x_size, y_size) ||
		!IsValidArmConfiguration(goalPos, numOfDOFs, map, x_size, y_size))
	{
		throw runtime_error("Invalid start or goal configuration!\n");
	}

	///////////////////////////////////////
	//// Feel free to modify anything below. Be careful modifying anything above.

	double **plan = NULL;
	int planlength = 0;

	// Call the corresponding planner function
	if (whichPlanner == PRM)
	{
		plannerPRM(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
	}
	else if (whichPlanner == RRT)
	{
		plannerRRT(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
	}
	else if (whichPlanner == RRTCONNECT)
	{
		plannerRRTConnect(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
	}
	else if (whichPlanner == RRTSTAR)
	{
		plannerRRTStar(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
	}
	else
	{
		planner(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
	}

	//// Feel free to modify anything above.
	//// If you modify something below, please change it back afterwards as my
	//// grading script will not work and you will recieve a 0.
	///////////////////////////////////////

	// Your solution's path should start with startPos and end with goalPos
	cout << "planlengthhtt " << planlength << endl;
	// int a=planlength;
	// 	for (int j=0;j<a;j++)
	// 	{
	// 		cout<<j<<endl;
	// 		for (int k=0;k<5;k++){
	// 			cout<<(plan)[j][k]<<" ";
	// 		}
	// 		cout<<"jjjjjjjjjjjjj = "<<j<<endl;
	// 	}

	if (!equalDoubleArrays(plan[0], startPos, numOfDOFs) ||
		!equalDoubleArrays(plan[planlength - 1], goalPos, numOfDOFs))
	{
		throw std::runtime_error("Start or goal position not matching");
	}

	/** Saves the solution to output file
	 * Do not modify the output log file output format as it is required for visualization
	 * and for grading.
	 */
	std::ofstream m_log_fstream;
	m_log_fstream.open(outputFile, std::ios::trunc); // Creates new or replaces existing file
	if (!m_log_fstream.is_open())
	{
		throw std::runtime_error("Cannot open file");
	}
	m_log_fstream << argv[1] << endl; // Write out map name first
	/// Then write out all the joint angles in the plan sequentially
	for (int i = 0; i < planlength; ++i)
	{
		for (int k = 0; k < numOfDOFs; ++k)
		{
			m_log_fstream << plan[i][k] << ",";
		}
		m_log_fstream << endl;
	}
}
