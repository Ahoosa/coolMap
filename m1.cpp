
#include "draw.h"
#include "m1.h"
#include "StreetsDatabaseAPI.h"
#include "OSMDatabaseAPI.h"
#include <iostream>
#include <cctype>
#include <algorithm>
#include <iterator>
#include <vector>
#include <numeric>
#include <map>
#include <set>
#include <cmath>
#include <unordered_map>
#include <math.h>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/range/algorithm_ext.hpp>


using namespace std;

//helper function
std::string changeName(std::string map_path);

//Global variables
int METERS_PER_KILOMETER = 1000;
int SECONDS_PER_HOUR = 3600;
double maxCitySpeed=0;

std::multimap<std::string, int> streetNameandID; //StreetNames and their relative IDs
std::unordered_map<OSMID, LatLon > nodeIDandLATLON; //NodeOSMID->LatLon
std::unordered_map<OSMID, std::vector<OSMID>> wayIDandnodeID; //wayOSMID -> vector of nodeOSMIDs
std::vector<std::vector<int>> mapSegments; //streetID->segID
std::vector<std::set<IntersectionIndex>> mapIntersections; //StreetID->Intersections
std::vector<std::vector<int>> intersectionSegments; //intersectionID->streetSegments
std::vector<std::pair<int, double>> mapSegmentsTravel; //Maps segments and their relative travel speeds

//Replacing .streets to .osm in map name

std::string changeName(std::string map_path) {
    std::string newName;
    std::size_t f = map_path.find(".streets");
    if (f == string::npos) {
        return " ";
    } else {
        return newName = map_path.replace(map_path.find(".streets"), sizeof (".streets") - 1, ".osm");
    }
}

bool load_map(std::string map_path) {

    std::string newName = changeName(map_path);
    bool load_successful = loadStreetsDatabaseBIN(map_path);
    bool OSMsuccess = loadOSMDatabaseBIN(newName);
     maxCitySpeed=0;

    //if unable to load maps:
    if (load_successful == false || OSMsuccess == false) {
        return false;
    } else {

        mapSegments.resize(getNumStreetSegments());
        mapIntersections.resize(getNumIntersections());
        intersectionSegments.resize(getNumIntersections());
        mapSegmentsTravel.resize(getNumStreetSegments());

        //filling vector of map Segments and relative travel time
        //Don't need to resize before starting loop because push_back() will do that for us
        for (StreetSegmentIndex i = 0; i < getNumStreetSegments(); i++) {
            InfoStreetSegment segment_info = getInfoStreetSegment(i);
            double travel_time = find_street_segment_length(i) / (segment_info.speedLimit * METERS_PER_KILOMETER / SECONDS_PER_HOUR);
            std::pair<int, double> pair1(i, travel_time);
            mapSegmentsTravel.push_back(pair1);
            if(segment_info.speedLimit>maxCitySpeed){
                maxCitySpeed=segment_info.speedLimit;
            }
        }

        //to find related data using the streetID (in order intersections and segments)
        for (StreetSegmentIndex i = 0; i < getNumStreetSegments(); ++i) {
            //Filling the mapSegment vector of vectors
            InfoStreetSegment segInfo = getInfoStreetSegment(i);
            mapSegments[segInfo.streetID].push_back(i);
            // Filling mapIntersection vector of vectors 
            IntersectionIndex fromInt = segInfo.from;
            IntersectionIndex toInt = segInfo.to;
            mapIntersections[segInfo.streetID].insert(fromInt);
            mapIntersections[segInfo.streetID].insert(toInt);

        }


        //vector of vectors to find segments of an intersection
        for (int inter = 0; inter < getNumIntersections(); ++inter) {
            for (int i = 0; i < getIntersectionStreetSegmentCount(inter); ++i) {
                int ss_id = getIntersectionStreetSegment(inter, i);
                intersectionSegments[inter].push_back(ss_id);
            }
        }

        if (streetNameandID.size() == 0) {
            for (int i = 0; i < getNumStreets(); ++i) {
                std::string street_name = getStreetName(i);
                //remove whitespace from street name and turn to lower case
                boost::algorithm::erase_all(street_name, " ");
                boost::algorithm::to_lower(street_name);
                //adding pairs of name and streetIDs to multimap
                std::pair<std::string, int> NameandID(street_name, i);
                streetNameandID.insert(NameandID);
            }
        }



        //unordered map that maps OSMway ids to the vector of the nodes that form the way
        for (int i = 0; i < getNumberOfWays(); ++i) {
            const OSMWay* way_osm = getWayByIndex(i);
            OSMID wayID = way_osm->id();
            //vector of ids of the nodes that make the way
            const std::vector<OSMID>& wayNode = getWayMembers(way_osm);
            std::pair<OSMID, std::vector < OSMID >> pairOfIDandNode(wayID, wayNode);
            wayIDandnodeID.insert(pairOfIDandNode);
        }

        //unordered map to map LATLON coordinates to nodeIDs  
        for (int i = 0; i < getNumberOfNodes(); ++i) {
            const OSMNode* node = getNodeByIndex(i);
            OSMID nodeID = node->id();
            LatLon coord = getNodeCoords(node);
            std::pair<OSMID, LatLon > pairOfIDandNode(nodeID, coord);
            nodeIDandLATLON.insert(pairOfIDandNode);
        }
        return load_successful;

    }
}



//Cleaning-up map related data structures 

void close_map() {
    mapSegments.clear();
    mapIntersections.clear();
    intersectionSegments.clear();
    streetNameandID.clear();
    wayIDandnodeID.clear();
    nodeIDandLATLON.clear();
    mapSegmentsTravel.clear();
    closeStreetDatabase();
    closeOSMDatabase();

}

//returns the distance between two coordinates in metres

double find_distance_between_two_points(std::pair<LatLon, LatLon> points) {
    double x1, x2, y1, y2;
    double lat_avg = (double) (points.first.lat() + points.second.lat()) / 2;

    x1 = DEGREE_TO_RADIAN * points.first.lon() * cos(DEGREE_TO_RADIAN * lat_avg);
    y1 = DEGREE_TO_RADIAN * points.first.lat();
    x2 = DEGREE_TO_RADIAN * points.second.lon() * cos(DEGREE_TO_RADIAN * lat_avg);
    y2 = DEGREE_TO_RADIAN * points.second.lat();

    double distance_between_points = EARTH_RADIUS_METERS * sqrt((pow(y2 - y1, 2))+(pow(x2 - x1, 2)));
    return distance_between_points;
}

//Returns the length of the given street segment in meters

double find_street_segment_length(int street_segment_id) {
    InfoStreetSegment segment_info = getInfoStreetSegment(street_segment_id);

    LatLon from_intersection_coordinate = getIntersectionPosition(segment_info.from);
    LatLon to_intersection_coordinate = getIntersectionPosition(segment_info.to);

    // If there are *no* curve points, return the distance between the starting point and ending point        
    if (segment_info.curvePointCount == 0) {
        std::pair <LatLon, LatLon> points1(from_intersection_coordinate, to_intersection_coordinate);
        return find_distance_between_two_points(points1);
    }// If there *are* curve points, traverse through each curve point to add up the length of the entire path
    else {
        double distance_curve_points = 0.0;
        // Add up the distance between each curve point
        for (int i = 0; i < segment_info.curvePointCount - 1; i++) {
            std::pair<LatLon, LatLon> points(getStreetSegmentCurvePoint(i, street_segment_id), getStreetSegmentCurvePoint(i + 1, street_segment_id));
            distance_curve_points += find_distance_between_two_points(points);
        }
        //need distance of two more segments, from starting intersection to first curve point and from last curve point to end intersection
        std::pair<LatLon, LatLon> points3(from_intersection_coordinate, getStreetSegmentCurvePoint(0, street_segment_id));
        std::pair<LatLon, LatLon> points4(getStreetSegmentCurvePoint(segment_info.curvePointCount - 1, street_segment_id), to_intersection_coordinate);
        // Sum total street segment length (distance form all curve points plus distance from first and last curve points to respective intersections)
        double street_segment_length = distance_curve_points + find_distance_between_two_points(points3) + find_distance_between_two_points(points4);

        return street_segment_length;
    }
}

//Returns the travel time to drive a street segment in seconds 
//(time = distance/speed_limit)

double find_street_segment_travel_time(int street_segment_id) {

    return find_street_segment_length(street_segment_id) / (getInfoStreetSegment(street_segment_id).speedLimit * METERS_PER_KILOMETER / SECONDS_PER_HOUR);;

}



//Returns the nearest intersection to the given position

int find_closest_intersection(LatLon my_position) {

    std::pair<double, IntersectionIndex> position_intersectionID;
    std::pair<double, IntersectionIndex> position_intersectionID2;
    std::map<double, int> IntersectionDistance;

    //Puts all intersection IDs and their distance with my_position in a map 
    for (StreetSegmentIndex i = 0; i < getNumStreetSegments(); ++i) {
        InfoStreetSegment infoSeg = getInfoStreetSegment(i);
        //Getting position of to and from 
        IntersectionIndex fromInter = infoSeg.from;
        IntersectionIndex toInter = infoSeg.to;
        LatLon intersectionPos1 = getIntersectionPosition(fromInter);
        LatLon interectionPos2 = getIntersectionPosition(toInter);

        //the second position of the pairs have the IDs  
        position_intersectionID.second = fromInter;
        position_intersectionID2.second = toInter;

        //finding relative distance between my position and each of from and to intersections 
        std::pair<LatLon, LatLon> findDistance(intersectionPos1, my_position);
        std::pair<LatLon, LatLon> findDistance2(interectionPos2, my_position);
        position_intersectionID.first = find_distance_between_two_points(findDistance);
        position_intersectionID2.first = find_distance_between_two_points(findDistance2);

        //The map includes the distance and the relative ID
        IntersectionDistance.insert(position_intersectionID);
        IntersectionDistance.insert(position_intersectionID2);
    }
    //first one would be the smallest distance 
    auto it = IntersectionDistance.begin();
    return (*it).second;
}

//Returns the street segments for the given intersection 

std::vector<int> find_street_segments_of_intersection(int intersection_id) {

    return intersectionSegments[intersection_id];

}

//Returns the street names at the given intersection (includes duplicate street 
//names in returned vector)

std::vector<std::string> find_street_names_of_intersection(int intersection_id) {

    std::vector<std::string> stNames;
    std::vector<int> streetSegs;
    int numSeg = getIntersectionStreetSegmentCount(intersection_id);
    // Store all of the intersection street segment ids in streetSegs vector
    for (int i = 0; i < numSeg; ++i) {
        int ss_id = getIntersectionStreetSegment(intersection_id, i);
        streetSegs.push_back(ss_id);
    }
    // Get the street segment info from the street segment id and 
    // store the street name string in the stNames vector
    for (auto& i : streetSegs) {
        InfoStreetSegment segment_information = getInfoStreetSegment(i);
        std::string segmentName = getStreetName(segment_information.streetID);
        stNames.push_back(segmentName);

    }
    return stNames;

}

//Returns true if you can get from intersection_ids.first to intersection_ids.second using a single 
//street segment (hint: check for 1-way streets too)
//corner case: an intersection is considered to be connected to itself

bool are_directly_connected(std::pair<int, int> intersection_ids) {
    // Store each intersection id in a separate vector
    std::vector<int> inter1Segments = intersectionSegments[intersection_ids.first];
    std::vector<int> inter2Segments = intersectionSegments[intersection_ids.second];
    // Sort intersection ids from each intersection in ascending order
    std::sort(inter1Segments.begin(), inter1Segments.end());
    std::sort(inter2Segments.begin(), inter2Segments.end());

    std::vector<int> commonSegment;
    std::set_intersection(inter1Segments.begin(),
            inter1Segments.end(),
            inter2Segments.begin(),
            inter2Segments.end(),
            std::back_inserter(commonSegment));
    if (intersection_ids.first == intersection_ids.second) {
        return true;
    } else if (commonSegment.size() == 0) {
        return false;
    } else {
        return true;
    }
}

//Returns all intersections reachable by traveling down one street segment 
//from given intersection (hint: you can't travel the wrong way on a 1-way street)
//the returned vector should NOT contain duplicate intersections

std::vector<int> find_adjacent_intersections(int intersection_id) {
    std::vector<int> adj;
    std::vector<int> adjStreetSegs = find_street_segments_of_intersection(intersection_id); //returns adj street seg IDs
    // Check each adjacent street segment and find intersection id of connecting intersection
    // if it can be reached by traveling from the current intersection
    for (auto i : adjStreetSegs) {
        InfoStreetSegment segment_information = getInfoStreetSegment(i);
        // For a one way street - only store intersection id if street segment in correct direction
        if (segment_information.oneWay) {
            if ((segment_information.from == intersection_id))
                adj.push_back(segment_information.to);
            // For a two-way street - push_back either from or to index
        } else if ((segment_information.from == intersection_id)) {
            adj.push_back(segment_information.to);
        } else {
            adj.push_back(segment_information.from);
        }
    }
    // Check for duplicate entries and erase as necessary
    std::sort(adj.begin(), adj.end());
    adj.erase(unique(adj.begin(), adj.end()), adj.end());
    return adj;
}


//Returns all street segments for the given street

std::vector<int> find_street_segments_of_street(int street_id) {
    // Store all map segments for the given street id and erase duplicates
    std::vector<int> stSegs = mapSegments[street_id];
    std::sort(stSegs.begin(), stSegs.end());
    stSegs.erase(unique(stSegs.begin(), stSegs.end()), stSegs.end());
    return stSegs;

}

//Returns all intersections along the a given street

std::vector<int> find_intersections_of_street(int street_id) {
    std::vector<int> stInts;
    // Ensure enough space for all map intersections
    stInts.reserve(mapIntersections[street_id].size());
    // Insert each map intersection into stInts vector
    for (auto val : mapIntersections[street_id]) {
        stInts.emplace_back(val);
    }
    return stInts;
}



//Return all intersection ids for two intersecting streets
//This function will typically return one intersection id.

std::vector<int> find_intersections_of_two_streets(std::pair<int, int> street_ids) {
    // Insert all intersections of each street into separate vectors
    std::vector<int> streetSegs1 = find_intersections_of_street(street_ids.first);
    std::vector<int> streetSegs2 = find_intersections_of_street(street_ids.second);

    std::vector<int> commonIntersection;
    // Find intersection ids that are common between both vectors 
    std::set_intersection(streetSegs1.begin(),
            streetSegs1.end(),
            streetSegs2.begin(),
            streetSegs2.end(),
            std::inserter(commonIntersection, commonIntersection.begin()));

    return commonIntersection;
}

//Returns all street ids corresponding to street names that start with the given prefix
//The function should be case-insensitive to the street prefix. You should ignore spaces.
//For example, both "bloor " and "BloOrst" are prefixes to "Bloor Street East".
//If no street names match the given prefix, this routine returns an empty (length 0) 
//vector.
//You can choose what to return if the street prefix passed in is an empty (length 0) 
//string, but your program must not crash if street_prefix is a length 0 string.
 
std::vector<int> find_street_ids_from_partial_street_name(std::string street_prefix) {
    std::vector<int> stIDs;

    // If prefix is an empty string, there are no matching street ids to return
    if (street_prefix.size() == 0) {
        std::vector<int> nothing;
        return nothing;  
    } else {
        // Adjust prefix string to accurately search by removing whitespace and replacing capital letters
        boost::algorithm::erase_all(street_prefix, " ");
        boost::algorithm::to_lower(street_prefix);
        //setting up the upper bound for the search
        std::string upPre = street_prefix;
        upPre.back() = street_prefix[street_prefix.length() - 1] + 1;

        // Search within bounds for matching street names, insert matches into stIDs
        for (auto itr = streetNameandID.lower_bound(street_prefix); itr != streetNameandID.upper_bound(upPre); ++itr) {

            stIDs.push_back(itr->second);
        }
        return stIDs;
    }
}

//Returns the area of the given closed feature in square meters
//Assume a non self-intersecting polygon (i.e. no holes)
//Return 0 if this feature is not a closed polygon.

double find_feature_area(int feature_id) {
    int feature_point_count = getFeaturePointCount(feature_id);
    LatLon first_feature_point = getFeaturePoint(0, feature_id);
    LatLon last_feature_point = getFeaturePoint(feature_point_count - 1, feature_id);

    double area = 0.0;
    std::vector<double> latVal;
    std::vector<std::pair<double, double>> coordinates;
    coordinates.resize(feature_point_count);

    // If the first and last points do not match, feature is not a closed polygon
    if (first_feature_point.lat() != last_feature_point.lat() || first_feature_point.lon() != last_feature_point.lon()) {
        return 0;
    }// It is a polygon
    else {
        // Insert all feature points into latVal vector
        for (int i = 0; i < feature_point_count; i++) {
            LatLon point = getFeaturePoint(i, feature_id);
            latVal.push_back(point.lat());
        }
        // Find average value of latVal vector
        double avg = std::accumulate(latVal.begin(), latVal.end(), 0.0) / latVal.size();
        double x1, y1;
        // Convert longitude and latitude points into radians
        // Insert into xVal vector
        for (int i = 0; i < feature_point_count; i++) {
            LatLon points = getFeaturePoint(i, feature_id);
            x1 = DEGREE_TO_RADIAN * points.lon() * cos(DEGREE_TO_RADIAN * avg);
            y1 = DEGREE_TO_RADIAN * points.lat();

            coordinates[i].first = x1;
            coordinates[i].second = y1;
        }
        int j = feature_point_count - 1;
        // Add up areas using feature points in xVal
        for (int i = 0; i < feature_point_count; i++) {
            std::pair<double, double> p0 = coordinates[i];
            std::pair<double, double> p1 = coordinates[j];
            area += (p0.first + p1.first)*(p1.second - p0.second);
            j = i;
        }
        // Area calculation using formula
        return abs(area * 0.5 * EARTH_RADIUS_METERS * EARTH_RADIUS_METERS);
    }


}

//Returns the length of the OSMWay that has the given OSMID, in meters.
//To implement this function you will have to  access the OSMDatabaseAPI.h 
//functions.

double find_way_length(OSMID way_id) {

    double length = 0.0;
    std::unordered_map<OSMID, std::vector < OSMID>>::iterator it = wayIDandnodeID.find(way_id);
    //If ID found get the vector of NodeIDs and match them with their Lat/Lon points that are in another map
    if (it != wayIDandnodeID.end()) {
        std::vector<OSMID> NodeIds = it->second;
        std::vector<LatLon> allPoints;

        for (auto x : NodeIds) {
            LatLon point = (nodeIDandLATLON.find(x))->second;
            allPoints.push_back(point);
        }
        //Once you have all the Lat/Lon points find the distance between all of them to get total length
        int j = 0;
        for (int i = 1; i < allPoints.size(); i++) {
            j = i - 1;
            length += find_distance_between_two_points(std::make_pair(allPoints[j], allPoints[i]));
        }
    }
    return length;
}
