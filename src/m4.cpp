using namespace std;


#include "m1.h"
#include "m2.h"
#include "m3.h"
#include "m4.h"
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


std::vector<int> currentInd;

int currentLoc[2], nextLoc[2];

int countDel[2];


int findFirstDepot(std::vector<DeliveryInfo> deliveries, std::vector<int> depots);

bool isLegal(std::vector<CourierSubpath> courier, const std::vector<int>& depots);
double heuristic3(int next, int end);
int bestStartingDepot1(std::vector<int> depots, std::vector<DeliveryInfo> deliveries);
int bestStartingDepot2(std::vector<int> depots, std::vector<DeliveryInfo> deliveries);
int secondBestStartingDepot2(std::vector<int> depots, std::vector<DeliveryInfo> deliveries);
int closestDepot(int currP, std::vector<int> depots);
double courierPathTravelTime(std::vector<CourierSubpath> courier, double turn_penalty);
std::vector<CourierSubpath> traveling_courier_firstD(
        const std::vector<DeliveryInfo>& deliveries,
        const std::vector<int>& depots,
        const float turn_penalty,
        const float truck_capacity,
        int firstD, int pathNumber);

// To choose a starting depot,
// will return the depot that has the shortest geometric distance to a pickup point
int bestStartingDepot1(std::vector<int> depots, std::vector<DeliveryInfo> deliveries) { //m*n 
    double smallestDistance = std::numeric_limits<double>::infinity();
    int bestDepot = depots[0];
    double distance = smallestDistance;
    for (int i = 0; i < depots.size(); ++i) {
        for (int j = 0; j < deliveries.size(); ++j) {
            distance = heuristic3(depots[i], deliveries[j].pickUp);
            if (distance < smallestDistance) bestDepot = depots[i];
        }
    }
    return bestDepot;
}
// To choose a starting depot,
// return the depot that has the most nearby pickup points (greedy)

//associate a score with each depot, return the one with the best score
// calculate the score by summing the geometric distances to each pickup point
// lowest score wins
int bestStartingDepot2(std::vector<int> depots, std::vector<DeliveryInfo> deliveries) {
    double bestDepotScore = std::numeric_limits<double>::infinity();
    //double depotScore = 0;
    int bestDepot = 0;
    for (int i = 0; i < depots.size(); ++i) {
        double depotScore = 0;
        for (int j = 0; j < deliveries.size(); ++j) {
            depotScore += heuristic3(depots[i], deliveries[j].pickUp);
        }
        if (depotScore < bestDepotScore) bestDepot = depots[i];
    }
    return bestDepot;
}
//return the second best starting depot, using same criteria as bestStartingDepot2
int secondBestStartingDepot2(std::vector<int> depots, std::vector<DeliveryInfo> deliveries) {
    double bestDepotScore = std::numeric_limits<double>::infinity();
    double secondBestScore = bestDepotScore;
    //double depotScore = 0;
    int bestDepot = 0, secondBest = 0;
    for (int i = 0; i < depots.size(); ++i) {
        double depotScore = 0;
        for (int j = 0; j < deliveries.size(); ++j) {
            depotScore += heuristic3(depots[i], deliveries[j].pickUp);
        }
        if (depotScore < bestDepotScore) {
            secondBest = bestDepot;
            bestDepot = depots[i];
        }
    }
    return secondBest;
    
}


double heuristic3(int next, int end) {

    return 3.6 * find_distance_between_two_points(std::make_pair(getIntersectionPosition(next), getIntersectionPosition(end))) / maxCitySpeed;
}

int closestDepot(int currP, std::vector<int> depots) {
    double distance = std::numeric_limits<double>::infinity();
    int depotP;
    int remove;

    for (int i = 0; i < depots.size(); ++i) {
        double estimate = heuristic3(currP, depots[i]);

        if (estimate <= distance) {
            distance = estimate;
            depotP = depots[i];
            remove = i;
        }

    }
    return depotP;
}

// return the closest point 
int closestPoint(std::vector<DeliveryInfo> deliveries, 
        int firstD, 
        std::vector<int> &visited, 
        std::vector<unsigned> &ind, 
        float truck_capacity, 
        std::vector<int> &currentDeliveries, 
        std::vector<int> &droppedOff, int pathNumber) {
    
    double distance = std::numeric_limits<double>::infinity();
    int firstP;
    int pick, drop, remove;

    std::vector<int> multiPickUp, multiDrop;

    // find total weight of all current things on truck
    float currentWeight = 0;
    for (int j = 0; j < currentDeliveries.size(); ++j) {
        currentWeight = currentWeight + deliveries[currentDeliveries[j]].itemWeight;
    }
    
    for (int i = 0; i < deliveries.size(); ++i) {

        double estimate = heuristic3(firstD, deliveries[i].pickUp);
        double estimate2 = heuristic3(firstD, deliveries[i].dropOff);
      
        if (estimate2 <= distance && visited[i] == 1 && droppedOff[i] == 0) {

            distance = estimate2;
            firstP = deliveries[i].dropOff;
            drop = i;
            pick = -1;
            remove = i;

            
        } else if (estimate <= distance && visited[i] == 0 && (currentWeight + deliveries[i].itemWeight) <= truck_capacity) {

            distance = estimate;
            firstP = deliveries[i].pickUp;
            pick = i;
            drop = -1;
            remove = i;
        }
    }


    if (currentLoc[pathNumber]  == 1) {
        for (auto index : currentInd) {
            ind.push_back(index);
        }
    }
    if (pick != -1) {
        multiPickUp.push_back(pick);
        currentWeight+=deliveries[pick].itemWeight;
        for (int i = 0; i < deliveries.size(); ++i) {
            if (i == pick) continue;
            else if (visited[i] == 0 && deliveries[pick].pickUp == deliveries[i].pickUp && (currentWeight + deliveries[i].itemWeight) <= truck_capacity) {
                multiPickUp.push_back(i);
                currentWeight+=deliveries[i].itemWeight;
            }
        }
        nextLoc[pathNumber] = 1;
        currentInd.clear();
        for (auto index : multiPickUp) {
            currentDeliveries.push_back(index);
            visited[index] = 1;
            currentInd.push_back(index);
        }
        

    } else if (drop != -1) {
        multiDrop.push_back(drop);

        for (int i = 0; i < deliveries.size(); ++i) {
            if (i == drop) continue;
            else if (visited[i] == 1 && droppedOff[i] == 0 && deliveries[drop].dropOff == deliveries[i].dropOff) {
                multiDrop.push_back(i);
            }
        }
        nextLoc[pathNumber] = 0;


        for (int index = 0; index < currentDeliveries.size(); ++index) {
            for (auto dropIndex : multiDrop) {
                if (currentDeliveries[index] == dropIndex) {
                    droppedOff[dropIndex] = 1;
                    countDel[pathNumber]--;
                    currentDeliveries.erase(currentDeliveries.begin() + index);
                }
            }
        }

    }

    multiPickUp.clear();
    multiDrop.clear();

   

    return firstP;
}

CourierSubpath subPathMaker(int firstD, int nextP, float turn_penalty, std::vector<unsigned>& ind, int pathNumber) {
    CourierSubpath info;
    std::vector<unsigned> empty;
    info.start_intersection = firstD;
    info.end_intersection = nextP;
    info.subpath = find_path_between_intersections(firstD, nextP, turn_penalty);
    int next = nextLoc[pathNumber];
    int curr = currentLoc[pathNumber];
    //-1 for depot 1 for pickup and 0 for dropoff
    info.pickUp_indices = ind;
    ind.clear(); //global variable so clear it so that it can be used for the next time

    return info;
}
double courierPathTravelTime(std::vector<CourierSubpath> courier, double turn_penalty) {
    // double compute_path_travel_time(const std::vector<StreetSegmentIndex>& path, const double turn_penalty)
    // get a vector of StreetSegmentIndex 's from the given courier vector of CourierSubpath's
    // each CourierSubpath has std::vector<int> subpath which is the street seg ids
    // for all courierSubpaths in courier, sum compute_travel_time(CourierSubpath[i].subpath, turn_penalty)
    double travelTime = 0;
    for (auto& i : courier) {
        travelTime += compute_path_travel_time(i.subpath, turn_penalty);
    }return travelTime;
}    
//start at best and second best depot, choose shortest path
std::vector<CourierSubpath> traveling_courier(
        const std::vector<DeliveryInfo>& deliveries,
        const std::vector<int>& depots,
        const float turn_penalty,
        const float truck_capacity) {
    
    int firstD1 = bestStartingDepot2(depots, deliveries);
    int firstD2 = secondBestStartingDepot2(depots, deliveries);
    //int firstD2 = bestStartingDepot1(depots, deliveries);
    
    std::vector<CourierSubpath> courier1 = traveling_courier_firstD(deliveries, depots, turn_penalty, truck_capacity, firstD1, 0);
    std::vector<CourierSubpath> courier2 = traveling_courier_firstD(deliveries, depots, turn_penalty, truck_capacity, firstD2, 1);
    
    //if travel time for courier 1 < that of courier 2, return courier 1, else return 2    
    double time1 = courierPathTravelTime(courier1, turn_penalty);
    double time2 = courierPathTravelTime(courier2, turn_penalty);
    
    if (isLegal(courier1, depots)) {
        if (isLegal(courier2, depots) && (time2 < time1)) return courier2;
        else return courier1;
    }
    std::vector<CourierSubpath> dumbCourier;
    return dumbCourier;
    //return courier1;
}
// returns vector of courierSubpaths
std::vector<CourierSubpath> traveling_courier_firstD(
        const std::vector<DeliveryInfo>& deliveries,
        const std::vector<int>& depots,
        const float turn_penalty,
        const float truck_capacity,
        int firstD, int pathNumber) {
    
    currentLoc[pathNumber] = -1;
    nextLoc[pathNumber] = 1;

    // vector to store ints of visited deliveries and droppedOff deliveries (delivery indexes?)
    std::vector<int> visited;
    std::vector<int> droppedOff;
    // 
    droppedOff.resize(deliveries.size());
    visited.resize(deliveries.size());
    for (int i = 0; i < deliveries.size(); ++i) {
        visited[i] = 0;
        droppedOff[i] = 0;
    }
    
    std::vector<int> currentDeliveries;
    std::vector<unsigned> ind;
    std::vector<CourierSubpath> courier;

//    int firstD = depots[findFirstDepot(deliveries, depots)];
//    int nextP = closestPoint(deliveries, firstD, visited, ind, truck_capacity, currentDeliveries, droppedOff);
//    CourierSubpath info = subPathMaker(firstD, nextP, turn_penalty, ind);
    countDel[pathNumber] = deliveries.size();
    // first dropOff point (int) is the first depot in the depots vector ////////////////////////////////////////
    
    
    //int firstD = depots[0];
    //int bestStartingDepot1(std::vector<int> depots, std::vector<DeliveryInfo> deliveries) 
    //int firstD = bestStartingDepot2(depots, deliveries); //need to find the depot closest to the most number of pickup/drop off points
    
    // next pickUp is the return of the closestPoint function
    int nextP = closestPoint(deliveries, firstD, visited, ind, truck_capacity, currentDeliveries, droppedOff, pathNumber);
    // make a subpath
    CourierSubpath info = subPathMaker(firstD, nextP, turn_penalty, ind, pathNumber);
    // update current position
    int currentPos = nextP;
    // add the subpath into courier which is the path to return
    courier.push_back(info);
    // update current location
    currentLoc[pathNumber] = nextLoc[pathNumber];
    
    while (countDel[pathNumber] > 0) {
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        nextP = closestPoint(deliveries, currentPos, visited, ind, truck_capacity, currentDeliveries, droppedOff, pathNumber);
        info = subPathMaker(currentPos, nextP, turn_penalty, ind, pathNumber);
        currentPos = nextP;
        courier.push_back(info);
        currentLoc[pathNumber] = nextLoc[pathNumber];
    }

    int depotP=closestDepot(currentLoc[pathNumber],depots);
    
    info = subPathMaker(currentPos, depotP, turn_penalty, ind, pathNumber);

    courier.push_back(info);

    visited.clear();
    ind.clear();
    currentDeliveries.clear();

    return courier;
}
// returns true if courier path begins and ends at a depot
bool isLegal(std::vector<CourierSubpath> courier, const std::vector<int>& depots){
    bool legalStart = false, legalEnd = false;
    
    if (std::find(depots.begin(), depots.end(), courier[0].start_intersection) != depots.end())
        legalStart = true;
    if (std::find(depots.begin(), depots.end(), courier[courier.size()-1].end_intersection) != depots.end())
        legalEnd = true;
    
    return (legalStart && legalEnd);
}

//call at the start, find depot & pickup location that are closest to each other
int findFirstDepot(std::vector<DeliveryInfo> deliveries, std::vector<int> depots){
    int bestDepot;
    double distance = std::numeric_limits<double>::infinity();

    for(int i = 0; i < depots.size(); ++i){
        for(int j = 0; j < deliveries.size(); ++j){
            if(heuristic3(depots[i], deliveries[j].pickUp) < distance){
                distance = heuristic3(depots[i], deliveries[j].pickUp);
                bestDepot = i;
            }
        }
    }
    return bestDepot;
}

std::vector<int> findDeliveryOrder(const std::vector<DeliveryInfo>& deliveries,
        const std::vector<int>& depots,
        const float turn_penalty,
        const float truck_capacity){
}
