

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   draw.h
 * Author: saeifara
 *
 * Created on February 23, 2020, 12:35 PM
 */

#ifndef DRAW_H
#define DRAW_H
//#pragma once

#include "m2.h"
#include "m1.h"
#include "StreetsDatabaseAPI.h"
#include "OSMDatabaseAPI.h"
#include "ezgl/application.hpp"
#include "ezgl/graphics.hpp"
#include "ezgl/rectangle.hpp"
#include <set>

extern double maxCitySpeed;

enum HighWayType{
    motorway=0,
    trunk,
    primary,
    secondary,
    tertiary,
    unclassified,
    residential
};

//bool sortbysec(const std::pair<int,double> &a, const std::pair<int,double> &b); 

struct feature_data{
    bool closed_feature;
    std::string name;
    FeatureType type;
    std::vector<ezgl::point2d> featurePoints;
};

struct intersection_data{
    LatLon position;
    std::string name;
    bool highlight = false;
    bool highlight2 = false;
};

struct POI_data{
     LatLon position;
     std::string name;
     std::string type;
    // OSMID osmID;
     bool isAmenity=false;
     bool isResto=false;
     bool isBar=false;
     bool isPub=false;
     bool isFastFood=false;
     bool isCafe=false;
     bool isIce=false;
     

};
struct street_data{
    std::vector<int> streetSegments;  
    std::string name;
    std::string type;
};




std::string getSegmentType(int segIdx);
std::pair<int,ezgl::color> segmentSpecification(std::string segmentType);
void setColourandWidth(ezgl::renderer *g,std::string segType);
void drawStreet(ezgl::renderer *g,std::vector<int> segments);

void showSegNames(ezgl::application *app,  std::vector<int> segments,ezgl::point2d hoverP,bool &startDisplay);

ezgl::point2d getCentreOfSeg(LatLon fromPoint, LatLon toPoint);
float getSegRotation(LatLon fromPoint, LatLon toPoint);

std::string replaceAND(std::string intersectionName);
std::pair<std::string,std::string> getNameInfo(std::vector<intersection_data> &intersections);
std::string getInfo(std::vector<intersection_data> &intersections,std::vector<POI_data> &pointofinterest);
std::pair<std::pair<double,double>,std::pair<double,double>> getMaxMin();
double latAvg(double min_lat,double max_lat);
double x_from_lon(double lon);
double y_from_lat(double lat);
double lon_from_x(double x);
double lat_from_y(double y);
std::string find_city(std::string selectedCity);
void getPOIData();
void POIdraw(ezgl::renderer *g,std::vector<POI_data> &pointofinterest,double worldRatio,int &countPOI);

LatLon ll_from_p2d(ezgl::point2d p2d);
ezgl::point2d p2d_from_ll(LatLon ll);
void getOSMmap();
void clearData(std::vector<intersection_data> &intersections,std:: vector<feature_data> &features,std::vector<street_data> &streets,std::vector<POI_data> &pointofinterest);
void fillVec(std::vector<intersection_data> &intersections,std:: vector<feature_data> &features,std::vector<street_data> &streets,std::vector<POI_data> &pointofinterest );
void clearHighlight(std::vector<intersection_data> &intersections, std::vector<POI_data> &pointofinterest);


void poiGetName(ezgl::renderer *g, std::string POIname, std::string POItype, double worldRatio, ezgl::point2d textPoint);
void draw_arrow(ezgl::renderer* g, std::vector<int> segments);
void featureGetName(ezgl::renderer *g, std::string featureName, FeatureType featureType, double worldRatio, std::vector<ezgl::point2d> points,int fID);
double determineCount(double worldRatio);
void determinePOI(std::vector<POI_data> &pointofinterest);
ezgl::color featureSpecification(int featureType);


#endif /* DRAW_H */
