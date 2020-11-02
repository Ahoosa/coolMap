/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
using namespace std;


#include "m1.h"
#include "m2.h"
#include "m3.h"

#include <gdk/gdk.h>
#include "draw.h"
#include "camera.hpp"
#include "canvas.hpp"
#include "point.hpp"
#include "callback.hpp"
#include "StreetsDatabaseAPI.h"
#include "OSMDatabaseAPI.h"
#include "ezgl/application.hpp"
#include "ezgl/graphics.hpp"
#include "ezgl/rectangle.hpp"
#include <math.h>
#include <sstream>
#include <string>
#include <boost/algorithm/string.hpp>
#include <set>
#include <vector>
#include <algorithm>
#include <stdlib.h>

#include <boost/heap/priority_queue.hpp>
#include <forward_list>
#include <iterator>
#include <queue>
#include <list>


#define No_EDGE -1
#define euclidean(p1,p2) std::hypotf(p1.x-p2.x,p1.y-p2.y)
#define fixFactor 5000

double heuristic(int next, int end) {

    return 3.6 * find_distance_between_two_points(std::make_pair(getIntersectionPosition(next), getIntersectionPosition(end))) / maxCitySpeed;
}

double heuristic2(int next, int end, double walkSpeed) {

    return find_distance_between_two_points(std::make_pair(getIntersectionPosition(next), getIntersectionPosition(end))) / walkSpeed;
}

class Node {
public:
    int id; //intersectionID
    int reachingEdge; //previous id used to get to this intersection
    double bestTime; //best time so far

    Node(IntersectionIndex idxID, StreetSegmentIndex prevEdge) {
        id = idxID;
        reachingEdge = prevEdge;
        bestTime = std::numeric_limits<double>::infinity();
    }

    ~Node() {

    }

};

struct WaveElem {
    Node *node; //intersection 
    int edgeID; //street segment id
    double travelTime;

    WaveElem(Node *n, int id, double time) {
        node = n;
        edgeID = id;
        travelTime = time;
    }
};

std::vector<Node> nodeInts;
void addNodes();
void clearNodes();

std::pair<int, std::list<StreetSegmentIndex>> bfsTraceBack2(int destID);
std::pair<int, std::vector<StreetSegmentIndex>> find_path_between_intersections2(
        const IntersectionIndex intersect_id_start,
        const IntersectionIndex intersect_id_end,
        const double turn_penalty, std::vector<int> possibleInter);
bool bfsMe2(int destID, double turnPenalty, std::vector<int> possibleInter);
double segmentWalkTime(int segmentID, double walking_speed);
std::pair<std::vector<StreetSegmentIndex>, std::vector<StreetSegmentIndex>> getAll(int start_intersection, int end_intersection, double turn_penalty, double walking_speed, double walking_time_limit);
bool bfsWalk(Node* source, int destID, double turnPenalty, double walking_speed);
std::vector<StreetSegmentIndex> WalkTraceBack(int destID);
std::vector<int> getBestWalkIntersection(int start_intersection, int end_intersection, double turn_penalty, double walking_speed, double walking_time_limit);

std::vector<StreetSegmentIndex> bfsTraceBack(int destID);
double costSoFar(int currentSegment, int nextSegment, double turnPenalty);
bool bfsMe(Node* source, int destID, double turnPenalty);
bool penaltyCheck(int prevSeg, int currSeg);
std::vector<Node*> intersectionIDandNode;
std::vector<Node*> walkIntersectionIDandNode;
typedef std::pair<double, WaveElem> Pair;
typedef std::pair<double, std::pair<int, std::vector<int>>> Pair2;

class TurnType {
    std::string Left;
    std::string Right;
};

class comparing {
public:

    double operator()(const Pair p1, const Pair p2) {
        return p1.first >= p2.first;
    }
};

class comparing2 {
public:

    double operator()(const Pair2 p1, const Pair2 p2) {
        return p1.first <= p2.first;
    }
};

class comparing3 {
public:

    double operator()(const Pair p1, const Pair p2) {
        return p1.first >= p2.first;
    }
};

void giveDirections(std::vector<StreetSegmentIndex> path);
void drawPath(ezgl::renderer *g, std::vector<IntersectionIndex> path, IntersectionIndex fromID, IntersectionIndex toID, bool walk);
void giveDirectionsWalking(std::vector<IntersectionIndex> walk, std::vector<IntersectionIndex> drive);
void drawPathWalking(ezgl::renderer *g, std::vector<IntersectionIndex> path);
std::vector<double> crossProduct(std::vector<double> currStreet, std::vector<double> nextStreet);
bool leftTurn(int currentStreet, int nextStreet);




// Returns the time required to travel along the path specified, in seconds.
// The path is given as a vector of street segment ids, and this function can
// assume the vector either forms a legal path or has size == 0.  The travel
// time is the sum of the length/speed-limit of each street segment, plus the
// given turn_penalty (in seconds) per turn implied by the path.  If there is
// no turn, then there is no penalty. Note that whenever the street id changes
// (e.g. going from Bloor Street West to Bloor Street East) we have a turn.

double compute_path_travel_time(const std::vector<StreetSegmentIndex>& path,
        const double turn_penalty) {
    //vector of street segments will either be a legal path or will have size 0
    //travel time is time for each segment,plus turn penalty
    int numTurns = 0;
    double travelTime = 0;
    if (path.size() == 0) {
        return 0;
    }
    //go through every segment and compute travel time for it
    //also check if it's on the same street as the next index in the vector and if not add a turn penalty
    for (int pathIndex = 0; pathIndex < path.size(); ++pathIndex) {
        double segTravelTime = find_street_segment_travel_time(path[pathIndex]);
        if (pathIndex < path.size() - 1) {
            if (getInfoStreetSegment(path[pathIndex]).streetID != getInfoStreetSegment(path[pathIndex + 1]).streetID) {
                numTurns = numTurns + 1;
            }
        }
        travelTime = travelTime + segTravelTime;
    }
    return travelTime + (double) numTurns*turn_penalty;
}


// Returns a path (route) between the start intersection and the end
// intersection, if one exists. This routine should return the shortest path
// between the given intersections, where the time penalty to turn right or
// left is given by turn_penalty (in seconds).  If no path exists, this routine
// returns an empty (size == 0) vector.  If more than one path exists, the path
// with the shortest travel time is returned. The path is returned as a vector
// of street segment ids; traversing these street segments, in the returned
// order, would take one from the start to the end intersection.
//3 tests to pass:
//route must be legal - sequence of street segments to form connected path, respect 1-way streets, return emptry vector if no path exists
//route should have minimal travel time between two intersections
//find route quickly (performance)

std::vector<StreetSegmentIndex> find_path_between_intersections(
        const IntersectionIndex intersect_id_start,
        const IntersectionIndex intersect_id_end,
        const double turn_penalty) {
    addNodes();


    std::vector<StreetSegmentIndex> bestPath;
    Node * source = new Node(intersect_id_start, -1);

    bool exist = bfsMe(source, intersect_id_end, turn_penalty);
    if (exist) {
        bestPath = bfsTraceBack(intersect_id_end);
        //        bestPath.resize(pathFound.size());
        //        std::copy(pathFound.begin(), pathFound.end(),bestPath.begin());
    } else {
        delete source;
        clearNodes();
        return bestPath;
    }

    for (int i = 0; i < intersectionIDandNode.size(); ++i) {
        intersectionIDandNode[i] = NULL;
    }
    intersectionIDandNode.clear();


    delete source;
    clearNodes();
    return bestPath;
}

bool bfsMe(Node* source, int destID, double turnPenalty) {
    
    // "waveFront" is a priority_queue, largest values at bottom, smallest at top, objects have type Pair
    // Pair = std::pair<double, WaveElem>
    priority_queue<Pair, std::vector<Pair>, comparing> waveFront; 
    
    //adjust capacity of vector to be able to contain all intersections
    intersectionIDandNode.reserve(getNumIntersections()); 
    
    // create visited vector of boolean values, all false to begin
    // keeps track of which intersections have been visited
    std::vector<bool> visited(getNumIntersections(), false);
    // create distance vector of double values, all inf to begin
    // stores travel time to each intersection from source
    std::vector<double> distance(getNumIntersections(), std::numeric_limits<double>::infinity());

    //distance from source node to source id is zero
    distance[source->id] = 0; 
    //    std::vector<Node*> temp(getNumIntersections(), NULL);

    // add source element to waveFront, 0 is travel time?
    // WaveElem contains: Node *node (intersection), int edgeID (street segment id), double travelTime
    waveFront.push(std::make_pair(0, WaveElem(source, -1, 0)));
    // it iterator points to first node in nodeInts vector (vector of nodes)
    std::vector<Node>::iterator it = nodeInts.begin();
    
    // while the waveFront is not empty
    while (!(waveFront.empty())) {
        
        // topOflist is the top element of the waveFront
        std::pair<double, WaveElem> topOflist = waveFront.top();
        // second element in topOflist pair is the waveElem - curr is the WaveElem at the top
        WaveElem curr = topOflist.second; //get travel time of this
        // *currNode is a pointer to current node
        Node *currNode = curr.node;
        // remove top element of waveFront (which is what *currNode is right now)
        waveFront.pop();

        // if the curr travelTime is LESS than or equal to the bestTime for that node,
        // and it has not yet been visited, then ...
        if (currNode->bestTime >= curr.travelTime && !visited[currNode->id]) {
            // update intersectionIDandNode vector (just a vector of Node ptrs) at the node id to currNode
            intersectionIDandNode[currNode->id] = currNode;
            // update distance vector with newly found travelTime
            distance[currNode->id] = curr.travelTime;
            // update bestTime to the newly found travel time, which should be less than or equal to the bestTime
            currNode->bestTime = curr.travelTime;
            // mark the node as visited in visited vector
            visited[currNode->id] = true;
            // update reachingEdge with edge taken
            currNode->reachingEdge = curr.edgeID;
            
            // store all street segs of that intersection (node) in outEdge vector
            std::vector<int> outEdge = find_street_segments_of_intersection(curr.node->id);
            
            // if destination reached, return true
            if (currNode->id == destID) {
                return true;
            }
            // iterate through all elements of outEdge
            for (auto& segmentID : outEdge) {
                // check for one way street validity
                InfoStreetSegment info = getInfoStreetSegment(segmentID);
                if ((info.oneWay && info.from != currNode->id) || segmentID == curr.edgeID) {
                    continue;
                }
                // add turn penalty if necessary
                double penaltyTime = 0;
                if (currNode->reachingEdge != -1 && penaltyCheck(currNode->reachingEdge, segmentID)) {
                    penaltyTime = turnPenalty;
                }
                // update waveFront, depending on direction of graph on segment
                if (!visited[info.to] && info.to != currNode->id) {
                    Node *toNode = &(*(it + info.to));
                    toNode->reachingEdge = segmentID;
                    waveFront.push(std::make_pair(heuristic(info.to, destID) + penaltyTime + currNode->bestTime + find_street_segment_travel_time(segmentID), 
                            WaveElem(toNode, segmentID, penaltyTime + currNode->bestTime + find_street_segment_travel_time(segmentID))));
                    //                    temp.push_back(toNode);
                } else if (!visited[info.from] && info.from != currNode->id) {
                    Node *fromNode = &(*(it + info.from));
                    fromNode->reachingEdge = segmentID;
                    waveFront.push(std::make_pair(heuristic(info.from, destID) + penaltyTime + currNode->bestTime + find_street_segment_travel_time(segmentID), WaveElem(fromNode, segmentID, penaltyTime + currNode->bestTime + find_street_segment_travel_time(segmentID))));
                    //                    temp.push_back(fromNode);
                }
            }
            outEdge.clear();

        }
    }


    return false;
}

bool penaltyCheck(int prevSeg, int currSeg) {
    InfoStreetSegment info1 = getInfoStreetSegment(prevSeg);
    InfoStreetSegment info2 = getInfoStreetSegment(currSeg);
    return info1.streetID != info2.streetID;
}

double costSoFar(int currentSegment, int nextSegment, double turnPenalty) {
    double val = find_street_segment_travel_time(nextSegment);
    InfoStreetSegment info1 = getInfoStreetSegment(currentSegment);
    InfoStreetSegment info2 = getInfoStreetSegment(nextSegment);

    if (info1.streetID == info2.streetID) return val;
    else {
        return val += turnPenalty;
    }
}


// Returns a path (route) between the start intersection and the end
// intersection, if one exists. This routine should return the shortest path
// between the given intersections, where the time penalty to turn right or
// left is given by turn_penalty (in seconds).  If no path exists, this routine
// returns an empty (size == 0) vector.  If more than one path exists, the path
// with the shortest travel time is returned. The path is returned as a vector
// of street segment ids; traversing these street segments, in the returned
// order, would take one from the start to the end intersection.
//3 tests to pass:
//route must be legal - sequence of street segments to form connected path, respect 1-way streets, return emptry vector if no path exists
//route should have minimal travel time between two intersections
//find route quickly (performance)

std::vector<StreetSegmentIndex> bfsTraceBack(int destID) {
    std::list<StreetSegmentIndex> path;

    Node *currNode = intersectionIDandNode[destID];
    auto prevEdge = currNode->reachingEdge;

    while (prevEdge != -1) {
        InfoStreetSegment info1 = getInfoStreetSegment(prevEdge);
        path.push_front(prevEdge);
        if (info1.from != currNode->id) {
            currNode = intersectionIDandNode[info1.from];
        } else {
            currNode = intersectionIDandNode[info1.to];
        }
        prevEdge = currNode->reachingEdge;
    }
    std::vector<StreetSegmentIndex> path2(path.size());

    std::copy(path.begin(), path.end(), path2.begin());
    return path2;
}


// This is an "uber pool"-like function. The idea is to minimize driving travel
// time by walking to a pick-up intersection (within walking_time_limit secs)
// from start_intersection while waiting for the car to arrive.  While walking,
// you can ignore speed limits of streets and just account for given
// walking_speed [m/sec]. However, pedestrians should also wait for traffic
// lights, so turn_penalty applies to whether you are driving or walking.
// Walking in the opposite direction of one-way streets is fine. Driving is
// NOT!  The routine returns a pair of vectors of street segment ids. The first
// vector is the walking path starting at start_intersection and ending at the
// pick-up intersection. The second vector is the driving path from pick-up
// intersection to end_interserction.  Note that a start_intersection can be a
// pick-up intersection. If this happens, the first vector should be empty
// (size = 0).  If there is no driving path found, both returned vectors should
// be empty (size = 0). 
// If the end_intersection turns out to be within the walking time limit, 
// you can choose what to return, but you should not crash. If your user 
// interface can call this function for that case, the UI should handle
// whatever you return in that case.

std::pair<std::vector<StreetSegmentIndex>, std::vector<StreetSegmentIndex>>
find_path_with_walk_to_pick_up(
        const IntersectionIndex start_intersection,
        const IntersectionIndex end_intersection,
        const double turn_penalty,
        const double walking_speed,
        const double walking_time_limit) {
    addNodes();
    return getAll(start_intersection, end_intersection, turn_penalty, walking_speed, walking_time_limit);
}

double segmentWalkTime(int segmentID, double walking_speed) {
    return find_street_segment_length(segmentID) / walking_speed;
}

std::pair<std::vector<StreetSegmentIndex>, std::vector<StreetSegmentIndex>> getAll(int start_intersection, int end_intersection, double turn_penalty, double walking_speed, double walking_time_limit) {
    std::vector<int> walkPath;
    double driveTime = std::numeric_limits<double>::infinity();

    if (walking_time_limit == 0) {
        return std::make_pair(walkPath, find_path_between_intersections(start_intersection, end_intersection, turn_penalty));
    }
    std::vector<int> closeIntersections = getBestWalkIntersection(start_intersection, end_intersection, turn_penalty, walking_speed, walking_time_limit);
    //std::priority_queue<Pair2, std::vector<Pair2>, comparing2> driveTimeandPath;
    std::pair<int, std::vector<int>> driveTimeandPath;

    if (closeIntersections.empty()) {
        return std::make_pair(walkPath, find_path_between_intersections(start_intersection, end_intersection, turn_penalty));
    }
    //    for (auto& walkInter : closeIntersections) {
    std::pair<int, std::vector<int>> result = find_path_between_intersections2(start_intersection, end_intersection, turn_penalty, closeIntersections);
    std::vector<int> drivePath = result.second;
    double pathTime = compute_path_travel_time(drivePath, turn_penalty);
    if (pathTime <= driveTime) {
        //              std::pair<int, std::vector<int>> walkInterandDrivePath(walkInter, drivePath);
        driveTime = pathTime;
        driveTimeandPath.first = result.first;
        driveTimeandPath.second = drivePath;
    }

    //    }

    //    Pair2 bestResult = driveTimeandPath.top();
    int lastWalkInt = driveTimeandPath.first;
    drivePath = driveTimeandPath.second;
    if (drivePath.empty()) {
        clearNodes();
        return std::make_pair(walkPath, drivePath);
    }

    Node * walkInt = new Node(start_intersection, -1);
    bool works = bfsWalk(walkInt, lastWalkInt, turn_penalty, walking_speed);
    if (works) {
        walkPath = WalkTraceBack(lastWalkInt);
    }

    delete walkInt;
    clearNodes();
    return std::make_pair(walkPath, drivePath);
}

std::vector<int> getBestWalkIntersection(int start_intersection, int end_intersection, double turn_penalty, double walking_speed, double walking_time_limit) {

    double timeSoFar = 0;
    int current = start_intersection;
    std::vector<int> closeIntersections;
    std::list<int> closeList;
    std::vector<bool> visited(getNumIntersections(), false);
    // closeIntersections.resize(getNumIntersections());
    Node * source = new Node(start_intersection, -1);
    closeList.push_back(current);
    while (!closeList.empty()) {
        current = closeList.front();
        closeList.pop_front();
        if (!visited[current]) {
            visited[current] = true;
            timeSoFar = 0;
            std::vector<int> segments = find_street_segments_of_intersection(current);
            for (auto& segID : segments) {
                InfoStreetSegment info = getInfoStreetSegment(segID);
                if (info.from != current) {
                    bool works = bfsWalk(source, info.from, turn_penalty, walking_speed);
                    if (works) {
                        std::vector<int> path = WalkTraceBack(info.from);
                        timeSoFar = compute_path_walking_time(path, walking_speed, turn_penalty);
                        if (timeSoFar <= walking_time_limit) {
                            closeIntersections.push_back(info.from);
                            closeList.push_front(info.from);
                            if (info.from == end_intersection) {
                                return closeIntersections;
                            }
                        }
                    }
                } else if (info.to != current) {
                    bool works = bfsWalk(source, info.to, turn_penalty, walking_speed);
                    if (works) {
                        std::vector<int> path = WalkTraceBack(info.to);
                        timeSoFar = compute_path_walking_time(path, walking_speed, turn_penalty);
                        if (timeSoFar <= walking_time_limit) {
                            closeIntersections.push_back(info.to);
                            closeList.push_front(info.to);
                            if (info.from == end_intersection) {
                                return closeIntersections;
                            }
                        }
                    }
                }
            }
        }
    }
    delete source;
    return closeIntersections;
}


//if intersection is within walking distance then check for walking path


// Returns the time required to "walk" along the path specified, in seconds.
// The path is given as a vector of street segment ids. The vector can be of
// size = 0, and in this case, it the function should return 0. The walking
// time is the sum of the length/<walking_speed> for each street segment, plus
// the given turn penalty, in seconds, per turn implied by the path. If there
// is no turn, then there is no penalty.  As mentioned above, going from Bloor
// Street West to Bloor street East is considered a turn

double compute_path_walking_time(const std::vector<StreetSegmentIndex>& path,
        const double walking_speed,
        const double turn_penalty) {
    int numTurns = 0;
    double walkingTime = 0;
    if (path.size() == 0) {
        return 0;
    }
    //go through vector and compute walking speed of every segment
    //if segment is on different street than next one in the vector, add 1 to numTurns
    for (int pathIndex = 0; pathIndex < path.size(); ++pathIndex) {
        double segWalkingTime = find_street_segment_length(path[pathIndex]) / walking_speed;
        if (pathIndex < path.size() - 1) {
            InfoStreetSegment currentSegInfo = getInfoStreetSegment(path[pathIndex]);
            InfoStreetSegment nextSegInfo = getInfoStreetSegment(path[pathIndex + 1]);
            if (currentSegInfo.streetID != nextSegInfo.streetID) {
                numTurns = numTurns + 1;
            }
        }
        walkingTime = walkingTime + segWalkingTime;
    }
    return walkingTime + (double) numTurns*turn_penalty;
}

bool bfsWalk(Node* source, int destID, double turnPenalty, double walking_speed) {
    std::priority_queue<Pair, std::vector<Pair>, comparing3> waveFront;

    walkIntersectionIDandNode.resize(getNumIntersections());

    std::vector<bool> visited(getNumIntersections(), false); //keeps track of which intersection we have checked
    std::vector<double> distance(getNumIntersections(), std::numeric_limits<double>::infinity()); //keeps traveltime of each intersection to the source

    distance[source->id] = 0; //distance from source to source is 0

    waveFront.push(std::make_pair(0, WaveElem(source, -1, 0)));

    std::vector<Node>::iterator it = nodeInts.begin();

    while (!(waveFront.empty())) {
        std::pair<double, WaveElem> topOflist = waveFront.top();
        WaveElem curr = topOflist.second; //get travel time of this
        //        double travelT = topOflist.first;
        Node *currNode = curr.node;
        waveFront.pop();


        if (!visited[currNode->id]) {
            walkIntersectionIDandNode[currNode->id] = currNode;
            distance[currNode->id] = curr.travelTime;
            currNode->bestTime = curr.travelTime;
            visited[currNode->id] = true;
            currNode->reachingEdge = curr.edgeID;
            std::vector<int> outEdge = find_street_segments_of_intersection(curr.node->id);

            if (currNode->id == destID) {
                return true;
            }

            for (auto& segmentID : outEdge) {

                InfoStreetSegment info = getInfoStreetSegment(segmentID);
                if (segmentID == curr.edgeID) {
                    continue;
                }
                double penaltyTime = 0;
                if (currNode->reachingEdge != -1 && penaltyCheck(currNode->reachingEdge, segmentID)) {
                    penaltyTime = turnPenalty;
                }
                if (!visited[info.to] && info.to != currNode->id) {

                    Node *toNode = &(*(it + info.to));

                    waveFront.push(std::make_pair(heuristic2(info.to, destID, walking_speed) + penaltyTime + currNode->bestTime + segmentWalkTime(segmentID, walking_speed), WaveElem(toNode, segmentID, penaltyTime + currNode->bestTime + segmentWalkTime(segmentID, walking_speed))));

                } else if (!visited[info.from] && info.from != currNode->id) {
                    Node *fromNode = &(*(it + info.from));

                    waveFront.push(std::make_pair(heuristic2(info.from, destID, walking_speed) + penaltyTime + currNode->bestTime + segmentWalkTime(segmentID, walking_speed), WaveElem(fromNode, segmentID, penaltyTime + currNode->bestTime + segmentWalkTime(segmentID, walking_speed))));

                }
            }
            outEdge.clear();

        }
    }

    return false;
}

std::vector<StreetSegmentIndex> WalkTraceBack(int destID) {
    std::list<StreetSegmentIndex> path;
    Node *currNode = walkIntersectionIDandNode[destID];
    auto prevEdge = currNode->reachingEdge;

    while (prevEdge != -1) {
        InfoStreetSegment info1 = getInfoStreetSegment(prevEdge);
        path.push_front(prevEdge);
        if (info1.from != currNode->id) {
            currNode = walkIntersectionIDandNode[info1.from];
        } else {
            currNode = walkIntersectionIDandNode[info1.to];
        }
        prevEdge = currNode->reachingEdge;
    }
    std::vector<int> vectPath;
    for (auto& c : path) {
        vectPath.push_back(c);
    }
    for (int i = 0; i < walkIntersectionIDandNode.size(); ++i) {
        walkIntersectionIDandNode[i] = NULL;
    }
    walkIntersectionIDandNode.clear();
    return vectPath;
}

std::pair<int, std::vector<StreetSegmentIndex>> find_path_between_intersections2(
        const IntersectionIndex intersect_id_start,
        const IntersectionIndex intersect_id_end,
        const double turn_penalty, std::vector<int> possibleInter) {

    std::vector<StreetSegmentIndex> bestPath;
    std::pair<int, std::list<int>> pathFound;

    bool exist = bfsMe2(intersect_id_end, turn_penalty, possibleInter);
    if (exist) {
        pathFound = bfsTraceBack2(intersect_id_end);
        for (auto &c : pathFound.second) {
            bestPath.push_back(c);
        }
    } else {
        return std::make_pair(intersect_id_start, bestPath);
    }
    for (int i = 0; i < intersectionIDandNode.size(); ++i) {
        intersectionIDandNode[i] = NULL;
    }
    intersectionIDandNode.clear();


    return std::make_pair(pathFound.first, bestPath);

}

bool bfsMe2(int destID, double turnPenalty, std::vector<int> possibleInter) {
    std::priority_queue<Pair, std::vector<Pair>, comparing> waveFront;
    intersectionIDandNode.reserve(getNumIntersections());

    std::vector<bool> visited(getNumIntersections(), false); //keeps track of which intersection we have checked

    //    std::vector<double> distance(getNumIntersections(), std::numeric_limits<double>::infinity()); //keeps traveltime of each intersection to the source
    Node *source = NULL;


    std::vector<Node>::iterator it = nodeInts.begin();

    for (auto& intersect : possibleInter) {
        source = &(*(it + intersect));
        source->reachingEdge = -1;
        waveFront.push(std::make_pair(0, WaveElem(source, -1, 0)));
    }

    while (!(waveFront.empty())) {
        std::pair<double, WaveElem> topOflist = waveFront.top();
        WaveElem curr = topOflist.second; //get travel time of this
        //        double travelT = topOflist.first;
        Node *currNode = curr.node;
        waveFront.pop();


        if (currNode->bestTime >= curr.travelTime && !visited[currNode->id]) {
            intersectionIDandNode[currNode->id] = currNode;
            //            distance[currNode->id] = curr.travelTime;
            currNode->bestTime = curr.travelTime;
            visited[currNode->id] = true;
            currNode->reachingEdge = curr.edgeID;
            std::vector<int> outEdge = find_street_segments_of_intersection(curr.node->id);

            if (currNode->id == destID) {
                return true;
            }

            for (auto& segmentID : outEdge) {

                InfoStreetSegment info = getInfoStreetSegment(segmentID);
                if ((info.oneWay && info.from != currNode->id) || segmentID == curr.edgeID) {
                    continue;
                }
                double penaltyTime = 0;
                if (currNode->reachingEdge != -1 && penaltyCheck(currNode->reachingEdge, segmentID)) {
                    penaltyTime = turnPenalty;
                }
                if (!visited[info.to] && info.to != currNode->id) {
                    Node *toNode = &(*(it + info.to));

                    waveFront.push(std::make_pair(heuristic(info.to, destID) + penaltyTime + currNode->bestTime + find_street_segment_travel_time(segmentID), WaveElem(toNode, segmentID, penaltyTime + currNode->bestTime + find_street_segment_travel_time(segmentID))));

                } else if (!visited[info.from] && info.from != currNode->id) {
                    Node *fromNode = &(*(it + info.from));

                    waveFront.push(std::make_pair(heuristic(info.from, destID) + penaltyTime + currNode->bestTime + find_street_segment_travel_time(segmentID), WaveElem(fromNode, segmentID, penaltyTime + currNode->bestTime + find_street_segment_travel_time(segmentID))));

                }
            }
            outEdge.clear();

        }
    }
    delete source;

    return false;
}

std::pair<int, std::list<StreetSegmentIndex>> bfsTraceBack2(int destID) {
    std::list<StreetSegmentIndex> path;
    Node *currNode = intersectionIDandNode[destID];
    auto prevEdge = currNode->reachingEdge;

    while (prevEdge != -1) {
        InfoStreetSegment info1 = getInfoStreetSegment(prevEdge);
        path.push_front(prevEdge);
        if (info1.from != currNode->id) {
            currNode = intersectionIDandNode[info1.from];
        } else {
            currNode = intersectionIDandNode[info1.to];
        }
        prevEdge = currNode->reachingEdge;
    }

    return std::make_pair(currNode->id, path);
}
void addNodes(){
    for(int i =0; i<getNumIntersections(); ++i){
        nodeInts.push_back(Node(i, 0));
    }
}
void clearNodes(){
    nodeInts.clear();
}

void drawPath(ezgl::renderer *g, std::vector<int> path, IntersectionIndex fromID, IntersectionIndex toID, bool walk) {

    //draw each segment on path in red. Also highlight start intersection in green and end intersection in end
    for (int pathIndex = 0; pathIndex < path.size(); ++pathIndex) {

        g->set_color(ezgl::LIGHT_SKY_BLUE);
        std::vector<int> segment;
        segment.push_back(path[pathIndex]);
        drawStreet(g, segment);

        double x_coord_from = x_from_lon(getIntersectionPosition(fromID).lon());
        double y_coord_from = y_from_lat(getIntersectionPosition(fromID).lat());
        double x_coord_to = x_from_lon(getIntersectionPosition(toID).lon());
        double y_coord_to = y_from_lat(getIntersectionPosition(toID).lat());
        g->set_color(ezgl::YELLOW);

        g->fill_arc({x_coord_to, y_coord_to}, 0.000005, 0, 360); 

        if (!walk) {
            g->set_color(ezgl::YELLOW);
            g->fill_arc({x_coord_from, y_coord_from}, 0.000005, 0, 360);
            g->set_color(ezgl::BLACK);
            g->draw_text({x_coord_from, y_coord_from}, "Starting point");
        }
        g->set_color(ezgl::BLACK);
        g->draw_text({x_coord_to, y_coord_to}, "Destination");
    }
}

//prints directions to terminal
void giveDirections(std::vector<StreetSegmentIndex> path) {
    bool North, South, East, West;
    North = false;
    South = false;
    East = false;
    West = false;
    //given a path of street segment indexes, give instructions to either move "N,S,E,W etc", left/right, or continue
    for (int pathIndex = 0; pathIndex < path.size(); ++pathIndex) {
        //for first instruction give a direction (N,S,E,W)
        LatLon fromIntersection = getIntersectionPosition(getInfoStreetSegment(path[pathIndex]).from);
        LatLon toIntersection = getIntersectionPosition(getInfoStreetSegment(path[pathIndex]).to);
        
        North = toIntersection.lat() > fromIntersection.lat();
        South = fromIntersection.lat() > toIntersection.lat();
        East = toIntersection.lon() > fromIntersection.lon();
        West = fromIntersection.lon() > toIntersection.lon();
        if(pathIndex==0){
            if(North && East){
                std::cout <<"Drive North-East on " <<getStreetName(getInfoStreetSegment(path[pathIndex]).streetID) <<"\n";
            }   
            else if(North && West){
                std::cout <<"Drive North-West on " <<getStreetName(getInfoStreetSegment(path[pathIndex]).streetID) <<"\n";
            }
            else if(North){
                std::cout <<"Drive North on " <<getStreetName(getInfoStreetSegment(path[pathIndex]).streetID) <<"\n";
            }
            else if(South && East){
                 std::cout <<"Drive South-East on " <<getStreetName(getInfoStreetSegment(path[pathIndex]).streetID) <<"\n";
            }
            else if(South && West){
                std::cout <<"Drive South-West on " <<getStreetName(getInfoStreetSegment(path[pathIndex]).streetID) <<"\n";
            }
            else if(South && !West && !East){
                std::cout <<"Drive South on " <<getStreetName(getInfoStreetSegment(path[pathIndex]).streetID) <<"\n";
            }
            else if(!North && !South && East){
                std::cout <<"Drive East on " <<getStreetName(getInfoStreetSegment(path[pathIndex]).streetID) <<"\n";
            }
            else{
                std::cout <<"Drive West on " <<getStreetName(getInfoStreetSegment(path[pathIndex]).streetID) <<"\n";
            }
        }
        else{
            if(getInfoStreetSegment(path[pathIndex]).streetID == getInfoStreetSegment(path[pathIndex-1]).streetID){
                if(pathIndex%path.size()/5 ==0)
                    std::cout <<"Continue straight on " <<getStreetName(getInfoStreetSegment(path[pathIndex]).streetID) <<"\n";     
            }
            else {
              
                if(leftTurn(getInfoStreetSegment(path[pathIndex-1]).streetID, getInfoStreetSegment(path[pathIndex]).streetID)){
                    std::cout<< "Turn left onto " <<getStreetName(getInfoStreetSegment(path[pathIndex]).streetID) <<"\n";
                }
                else{
                    std::cout<< "Turn right onto " <<getStreetName(getInfoStreetSegment(path[pathIndex]).streetID) <<"\n";
                } 
            }
        }
    }
    std::cout<< "You have arrived at your destination\n";
}


//give walking directions
void giveDirectionsWalking(std::vector<IntersectionIndex> walk, std::vector<IntersectionIndex> drive) {
    bool North, South, East, West;
    North = false; South = false; West = false; East = false;
    if (walk.size() == 0) {
        std::cout << "You do not need to walk, you are already at a pickup intersection\n";
        giveDirections(drive);
    } else if (drive.size() == 0) {
        std::cout << "Your destination is within walking distance\n";
    }
    for(int walkIndex = 0; walkIndex < walk.size(); ++walkIndex){
        LatLon fromIntersection = getIntersectionPosition(getInfoStreetSegment(walk[walkIndex]).from);
        LatLon toIntersection = getIntersectionPosition(getInfoStreetSegment(walk[walkIndex]).to);
        North = toIntersection.lat() > fromIntersection.lat();
        South = fromIntersection.lat() > toIntersection.lat();
        East = toIntersection.lon() > fromIntersection.lon();
        West = fromIntersection.lon() > toIntersection.lon();
        if(walkIndex == 0){
            if(North && East){
                std::cout <<"Walk North-East on " <<getStreetName(getInfoStreetSegment(walk[walkIndex]).streetID) <<"\n";
            }   
            else if(North && West){
                std::cout <<"Walk North-West on " <<getStreetName(getInfoStreetSegment(walk[walkIndex]).streetID) <<"\n";
            }
            else if(North){
                std::cout <<"Walk North on " <<getStreetName(getInfoStreetSegment(walk[walkIndex]).streetID) <<"\n";
            }
            else if(South && East){
                 std::cout <<"Walk South-East on " <<getStreetName(getInfoStreetSegment(walk[walkIndex]).streetID) <<"\n";
            }
            else if(South && West){
                std::cout <<"Walk South-West on " <<getStreetName(getInfoStreetSegment(walk[walkIndex]).streetID) <<"\n";
            }
            else if(South && !West && !East){
                std::cout <<"Walk South on " <<getStreetName(getInfoStreetSegment(walk[walkIndex]).streetID) <<"\n";
            }
            else if(!North && !South && East){
                std::cout <<"Walk East on " <<getStreetName(getInfoStreetSegment(walk[walkIndex]).streetID) <<"\n";
            }
            else{
                std::cout <<"Walk West on " <<getStreetName(getInfoStreetSegment(walk[walkIndex]).streetID) <<"\n";
            }
        }
        else {
            if(getInfoStreetSegment(walk[walkIndex]).streetID == getInfoStreetSegment(walk[walkIndex-1]).streetID){
                std::cout<<"Continue Straight on " <<getStreetName(getInfoStreetSegment(walk[walkIndex]).streetID) <<"\n";
            }
            else{
                if(leftTurn(getInfoStreetSegment(walk[walkIndex-1]).streetID, getInfoStreetSegment(walk[walkIndex]).streetID)){
                    std::cout<< "Turn left onto " <<getStreetName(getInfoStreetSegment(walk[walkIndex]).streetID) <<"\n";
                }
                else{
                    std::cout<< "Turn right onto " <<getStreetName(getInfoStreetSegment(walk[walkIndex]).streetID) <<"\n";
                }
            }
        }    
    }
    if(drive.size()!=0){
        std::cout<< "You have arrived at the pickup point\n";
        giveDirections(drive);
    }
    else{
        std::cout<< "You have arrived at your destination\n";
    }
    

}


//draw walking path
void drawPathWalking(ezgl::renderer *g, std::vector<IntersectionIndex> path) {
    for (int pathIndex = 0; pathIndex < path.size(); ++pathIndex) {
        //        g->set_line_dash(ezgl::asymmetric_5_3);
        g->set_color(ezgl::LIGHT_PINK);
        std::vector<int> segment;
        segment.push_back(path[pathIndex]);
        drawStreet(g, segment);
        //        g->set_line_dash(none);

        double x_coord, y_coord;

        if (pathIndex == 0 && path.size() > 1) {
            if (getIntersectionName(getInfoStreetSegment(path[pathIndex]).from) == getIntersectionName(getInfoStreetSegment(path[pathIndex + 1]).from) || getIntersectionName(getInfoStreetSegment(path[pathIndex]).from) == getIntersectionName(getInfoStreetSegment(path[pathIndex + 1]).to)) {
                //to intersection of first segment is starting point 
                x_coord = x_from_lon(getIntersectionPosition(getInfoStreetSegment(path[pathIndex]).to).lon());
                y_coord = y_from_lat(getIntersectionPosition(getInfoStreetSegment(path[pathIndex]).to).lat());
            }
            else if (getIntersectionName(getInfoStreetSegment(path[pathIndex]).to) == getIntersectionName(getInfoStreetSegment(path[pathIndex + 1]).to) || getIntersectionName(getInfoStreetSegment(path[pathIndex]).to) == getIntersectionName(getInfoStreetSegment(path[pathIndex + 1]).from)) {
                //from intersection of first segment is starting point
                x_coord = x_from_lon(getIntersectionPosition(getInfoStreetSegment(path[pathIndex]).from).lon());
                y_coord = y_from_lat(getIntersectionPosition(getInfoStreetSegment(path[pathIndex]).from).lat());
            }
            g->set_color(ezgl::YELLOW);
            g->fill_arc({x_coord, y_coord}, 0.000005, 0, 360);
            g->set_color(ezgl::BLACK);
            g->draw_text({x_coord, y_coord}, "Starting Point");
        }
        if (pathIndex == (path.size() - 2) && path.size() > 1) {
            if (getIntersectionName(getInfoStreetSegment(path[pathIndex]).from) == getIntersectionName(getInfoStreetSegment(path[pathIndex + 1]).from) || getIntersectionName(getInfoStreetSegment(path[pathIndex]).to) == getIntersectionName(getInfoStreetSegment(path[pathIndex + 1]).from)) {
                //to intersection of last segment is pickup point 
                x_coord = x_from_lon(getIntersectionPosition(getInfoStreetSegment(path[pathIndex + 1]).to).lon());
                y_coord = y_from_lat(getIntersectionPosition(getInfoStreetSegment(path[pathIndex + 1]).to).lat());
            }
            else if (getIntersectionName(getInfoStreetSegment(path[pathIndex]).to) == getIntersectionName(getInfoStreetSegment(path[pathIndex + 1]).to) || getIntersectionName(getInfoStreetSegment(path[pathIndex]).from) == getIntersectionName(getInfoStreetSegment(path[pathIndex + 1]).to)) {
                //from intersection of last segment is pickup point 
                x_coord = x_from_lon(getIntersectionPosition(getInfoStreetSegment(path[pathIndex + 1]).from).lon());
                y_coord = y_from_lat(getIntersectionPosition(getInfoStreetSegment(path[pathIndex + 1]).from).lat());
            }
             g->set_color(ezgl::RED);
             g->fill_arc({x_coord, y_coord}, 0.000005, 0, 360);
             g->set_color(ezgl::BLACK);
             g->draw_text({x_coord, y_coord}, "Pickup Point");

        }
    }
}

std::vector<double> crossProduct(std::vector<double> currStreet, std::vector<double> nextStreet){
    std::vector<double> product;
    product.push_back(currStreet[1]*nextStreet[2]-nextStreet[1]*currStreet[2]);
    product.push_back(-1*(currStreet[0]*nextStreet[2] - currStreet[2]*nextStreet[0]));
    product.push_back(currStreet[0]*nextStreet[1] - currStreet[1]*nextStreet[0]);
    return product;
}

bool leftTurn(int currentStreet, int nextStreet){
    std::vector<double> currStreet;
    std::vector<double> nexStreet;
    std::vector<IntersectionIndex> intersection = find_intersections_of_two_streets(make_pair(currentStreet, nextStreet));
   
    currStreet.push_back(getIntersectionPosition(getInfoStreetSegment(currentStreet).from).lon()-getIntersectionPosition(intersection[0]).lon());
    currStreet.push_back(getIntersectionPosition(getInfoStreetSegment(currentStreet).from).lat()-getIntersectionPosition(intersection[0]).lat());
    currStreet.push_back(0);
    nexStreet.push_back(getIntersectionPosition(intersection[0]).lon()-getIntersectionPosition(getInfoStreetSegment(nextStreet).to).lon());
    nexStreet.push_back(getIntersectionPosition(intersection[0]).lat()-getIntersectionPosition(getInfoStreetSegment(nextStreet).to).lat());
    nexStreet.push_back(0);
    std::vector<double> product = crossProduct(currStreet, nexStreet);
    return (product[2]>0);
}
