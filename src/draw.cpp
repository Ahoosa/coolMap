

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


#include "m2.h"
#include "m1.h"
#include "draw.h"
#include "point.hpp"
#include "StreetsDatabaseAPI.h"
#include "OSMDatabaseAPI.h"
#include "ezgl/application.hpp"
#include "ezgl/graphics.hpp"
#include "ezgl/rectangle.hpp"
#include <cmath>

double lat_avg;
std::map<OSMID, std::string> wayIDandValue;
std::vector<std::pair<std::string, ezgl::point2d>> sameName;
std::map<OSMID, std::string> streetNum;
std::map<OSMID, std::string> streetName;

//get the segment type from the street segment index
std::string getSegmentType(StreetSegmentIndex segIdx) {
    InfoStreetSegment segInfo = getInfoStreetSegment(segIdx);
    OSMID wayID = segInfo.wayOSMID;
    std::map<OSMID, std::string>::iterator it = wayIDandValue.find(wayID);
    return it->second;
}

void getOSMmap() {
    //    wayIDandValue.reserve( getNumberOfWays());
    wayIDandValue.clear();
    for (int wayIdx = 0; wayIdx < getNumberOfWays(); ++wayIdx) {
        const OSMWay* e = getWayByIndex(wayIdx);
        for (int i = 0; i < getTagCount(e); ++i) {
            std::string key, value;
            std::tie(key, value) = getTagPair(e, i);
            if (key == "highway") {
                OSMID wayID2 = e->id();
                wayIDandValue.insert(std::make_pair(wayID2, value));
            }
        }

    }
    streetName.clear();
    streetNum.clear();
    for (int nIdx = 0; nIdx < getNumberOfNodes(); ++nIdx) {
        const OSMNode* n = getNodeByIndex(nIdx);
        for (int i = 0; i < getTagCount(n); ++i) {
            std::string key, value, key2, value2;
            std::tie(key, value) = getTagPair(n, i);
            std::tie(key2, value2) = getTagPair(n, i);
             OSMID nID2 = n->id();
            if (key == "addr:street" || key2== "addr:housenumber" ) {
               
                streetName.insert(std::make_pair(nID2, value));

            }  
             if (key2== "addr:housenumber" ) {
            
                streetNum.insert(std::make_pair(nID2, value2));
            }
        }

    }
   
}
//return a pair that contains an int and a specific colour, depending on the type
// of segment
std::pair<int, ezgl::color> segmentSpecification(std::string segmentType) {
    if (segmentType == "motorway") {
        return std::make_pair(6, ezgl::color(255, 140, 100, 255));
    } else if (segmentType == "trunk") {
        return std::make_pair(6, ezgl::PLUM);
    } else if (segmentType == "primary") {
        return std::make_pair(5, ezgl::ORANGE);
    } else if (segmentType == "secondary") {
        return std::make_pair(6, ezgl::BISQUE);
    } else if (segmentType == "tertiary") {
        return std::make_pair(5, ezgl::BISQUE);
    } else if (segmentType == "unclassified") {
        return std::make_pair(5, ezgl::WHITE);
    } else if (segmentType == "residential") {
        return std::make_pair(5, ezgl::WHITE);
    }

}

//specifies the colour and length of lines
void setColourandWidth(ezgl::renderer *g, std::string segType) {

    std::pair<double, ezgl::color> WidthandColor = segmentSpecification(segType);
    g->set_color(WidthandColor.second);
    g->set_line_width(WidthandColor.first);

}
// Draw the streets on the map
void drawStreet(ezgl::renderer *g, std::vector<int> segments) {
    for (auto segID : segments) {
        //retrieve necessary data using the segment id, for all segment ids contained within segments
        InfoStreetSegment segInfo = getInfoStreetSegment(segID);
        int segID2 = segID;
        LatLon intFrom = getIntersectionPosition(segInfo.from);
        LatLon intTo = getIntersectionPosition(segInfo.to);
        //if the street is within the visible range of the map, draw a straight or curved line tor represent the street
        if (g->get_visible_world().contains(p2d_from_ll(intFrom)) || g->get_visible_world().contains(p2d_from_ll(intTo))) {
            if (segInfo.curvePointCount == 0) {
                float x1 = x_from_lon(intFrom.lon());
                float y1 = y_from_lat(intFrom.lat());
                float x2 = x_from_lon(intTo.lon());
                float y2 = y_from_lat(intTo.lat());

                g->draw_line({x1, y1},
                {
                    x2, y2 });

            } else {
                for (int i = 0; i < segInfo.curvePointCount - 1; i++) {
                    LatLon curveP = getStreetSegmentCurvePoint(i, segID);
                    LatLon curveP2 = getStreetSegmentCurvePoint(i + 1, segID);
                    float x1 = x_from_lon(curveP.lon());
                    float y1 = y_from_lat(curveP.lat());
                    float x2 = x_from_lon(curveP2.lon());
                    float y2 = y_from_lat(curveP2.lat());
                    g->draw_line({x1, y1},
                    {
                        x2, y2 });

                }
                LatLon curveFirst = getStreetSegmentCurvePoint(0, segID);
                LatLon curveLast = getStreetSegmentCurvePoint(segInfo.curvePointCount - 1, segID);
                float x1 = x_from_lon(curveFirst.lon());
                float y1 = y_from_lat(curveFirst.lat());
                float x2 = x_from_lon(intFrom.lon());
                float y2 = y_from_lat(intFrom.lat());
                g->draw_line({x1, y1},
                {
                    x2, y2 });

                float x3 = x_from_lon(curveLast.lon());
                float y3 = y_from_lat(curveLast.lat());
                float x4 = x_from_lon(intTo.lon());
                float y4 = y_from_lat(intTo.lat());
                g->draw_line({x3, y3},
                {
                    x4, y4 });

            }
        }

    }
}
//display the known names of street segments when the user hovers their mouse over
void showSegNames(ezgl::application* app, std::vector<int> segments, ezgl::point2d hoverP, bool &startDisplay) {
    ezgl::renderer *g = app->get_renderer();
    g->set_color(ezgl::BLACK);
    g->format_font("monospace", ezgl::font_slant::normal, ezgl::font_weight::normal, 10);
    g->set_font_size(10);


    if (sameName.size() != 4 && startDisplay == false) {

        for (auto segID : segments) {
            InfoStreetSegment segInfo = getInfoStreetSegment(segID);

          if (g->get_visible_world().contains(p2d_from_ll(getIntersectionPosition(segInfo.from))) || g->get_visible_world().contains(p2d_from_ll(getIntersectionPosition(segInfo.to)))){
                        
            float segment_rotation = getSegRotation(getIntersectionPosition(segInfo.from), getIntersectionPosition(segInfo.to));
            g->set_text_rotation(segment_rotation);
            std::string seg_name = getStreetName(segInfo.streetID);


            startDisplay = false;
            if (g->get_visible_world().contains(p2d_from_ll(getIntersectionPosition(segInfo.from))) || g->get_visible_world().contains(p2d_from_ll(getIntersectionPosition(segInfo.to)))) {
                
                if (segInfo.curvePointCount == 0) {
                    ezgl::point2d text_location = getCentreOfSeg(getIntersectionPosition(segInfo.from), getIntersectionPosition(segInfo.to));
                    float x1 = x_from_lon(getIntersectionPosition(segInfo.from).lon());
                    float y1 = y_from_lat(getIntersectionPosition(segInfo.from).lat());
                    float x2 = x_from_lon(getIntersectionPosition(segInfo.to).lon());
                    float y2 = y_from_lat(getIntersectionPosition(segInfo.to).lat());
                    ezgl::rectangle textBox({x1, y1},
                    {  x2, y2 });
                    if (textBox.contains(hoverP)) {
                        if ((getStreetName(segInfo.streetID) != "<unknown>")) { 
                            sameName.push_back(std::make_pair(seg_name, text_location));
                        }
                    }
                } else {
                    for (int i = 0; i < segInfo.curvePointCount - 1; i+=2) {  
                        LatLon curveP = getStreetSegmentCurvePoint(i, segID);
                        LatLon curveP2 = getStreetSegmentCurvePoint(i + 1, segID);
                        ezgl::point2d text_location = getCentreOfSeg(curveP, curveP2);
                        float x1 = x_from_lon(curveP.lon());
                        float y1 = y_from_lat(curveP.lat());
                        float x2 = x_from_lon(curveP2.lon());
                        float y2 = y_from_lat(curveP2.lat());
                        ezgl::rectangle textBox({x1, y1},
                        { x2, y2  });

                        if (textBox.contains(hoverP)) {
                            if ((getStreetName(segInfo.streetID) != "<unknown>")) {
                                sameName.push_back(std::make_pair(seg_name, text_location));   
                            }
                        }
                    }

                }

            }
        }
        } 
    }else {
        startDisplay = true;
        if (std::adjacent_find(sameName.begin(), sameName.end(), std::not_equal_to<>()) == sameName.end()) {
            g->draw_text(sameName[0].second, sameName[0].first);   
        }
        sameName.clear();
        startDisplay = false;
    }
    app->flush_drawing();
}

//returns max and min bounds for a city
std::pair<std::pair<double, double>, std::pair<double, double>> getMaxMin() {
    double max_lat, min_lat, max_lon, min_lon;
    for (int id = 0; id < getNumIntersections(); ++id) {
        double cLat = getIntersectionPosition(id).lat();
        double cLon = getIntersectionPosition(id).lon();

        if (id == 0) {
            max_lat = cLat;
            min_lat = cLat;
            max_lon = cLon;
            min_lon = cLon;
        } else {
            if (cLat > max_lat) {
                max_lat = cLat;
            } else if (cLat < min_lat) {
                min_lat = cLat;
            }
            if (cLon > max_lon) {
                max_lon = cLon;
            } else if (cLon < min_lon) {
                min_lon = cLon;
            }
        }
    }
    std::pair<double, double> maxVal = std::make_pair(max_lat, max_lon);
    std::pair<double, double> minVal = std::make_pair(min_lat, min_lon);
    return std::make_pair(maxVal, minVal);
}

//find the selected city
std::string find_city(std::string selectedCity) {
    std::vector<std::string> cities = {"toronto_canada",
        "beijing_china",
        "cairo_egypt",
        "cape-town_south-africa",
        "golden-horseshoe_canada",
        "tehran_iran",
        "tokyo_japan",
        "sydney_australia",
        "rio-de-janeiro_brazil",
        "london_england",
        "new-york_usa",
        "new-delhi_india",
        "moscow_russia",
        "interlaken_switzerland",
        "iceland",
        "hong-kong_china",
        "hamilton_canada"};

    for (int i = 0; i < cities.size(); ++i) {
        std::size_t found = cities[i].find(selectedCity);
        if (found != std::string::npos) {
            return "/cad2/ece297s/public/maps/" + cities[i] + ".streets.bin";
        }
    }
}

//sorts features
bool sortbysec(const std::pair<int,double> &a, const std::pair<int,double> &b) { 
    return (a.second > b.second); 
}

//fill all the containers with city data
void fillVec(std::vector<intersection_data> &intersections, std::vector<feature_data> &features, std::vector<street_data> &streets, std::vector<POI_data> &pointofinterest) {

    getOSMmap();

    intersections.resize(getNumIntersections());
    features.resize(getNumFeatures());
    streets.resize(getNumStreets());
    pointofinterest.resize(getNumPointsOfInterest());

    for (int id = 0; id < getNumIntersections(); ++id) {
        intersections[id].position = getIntersectionPosition(id);
        intersections[id].name = getIntersectionName(id);
    }

    std::vector<std::pair<int, double>> indexArea;
    indexArea.resize(getNumFeatures());
    for(int featureIndex = 0; featureIndex <getNumFeatures(); ++featureIndex){
        indexArea[featureIndex].first = featureIndex;
        indexArea[featureIndex].second = find_feature_area(featureIndex);
    }
    
    std::sort(indexArea.begin(), indexArea.end(), sortbysec); 
    
     for (int featureID = 0; featureID < getNumFeatures(); ++featureID) {
        //fill vector with features in order form biggest area to smallest
        int index = indexArea[featureID].first;
        features[featureID].name = getFeatureName(index);
        features[featureID].type = getFeatureType(index);
         if ((getFeaturePoint(0, index).lat() == getFeaturePoint((getFeaturePointCount(index) - 1), index).lat()) &&
             (getFeaturePoint(0, index).lon() == getFeaturePoint((getFeaturePointCount(index) - 1), index).lon())) {
        features[featureID].closed_feature = true;
        }
        else {
            features[featureID].closed_feature = false;
        }
        for(int point = 0; point<getFeaturePointCount(index); ++point){
            features[featureID].featurePoints.push_back({x_from_lon(getFeaturePoint(point, index).lon()), y_from_lat(getFeaturePoint(point, index).lat())});
        }
     }

    for (int stID = 0; stID < getNumStreets(); ++stID) {
        std::vector<int> segs = find_street_segments_of_street(stID);
        streets[stID].streetSegments.resize(segs.size());
        streets[stID].streetSegments = segs;

        //         streets[stID].type.resize(segs.size());
        streets[stID].name = getStreetName(stID);

        streets[stID].type = getSegmentType(segs[0]);

        segs.clear();
    }



    for (int numPOI = 0; numPOI < getNumPointsOfInterest(); ++numPOI) {
        pointofinterest[numPOI].position = getPointOfInterestPosition(numPOI);
        pointofinterest[numPOI].name = getPointOfInterestName(numPOI);
        pointofinterest[numPOI].type = getPointOfInterestType(numPOI);

    }





}

//clear the data containers
void clearData(std::vector<intersection_data> &intersections, std::vector<feature_data> &features, std::vector<street_data> &streets, std::vector<POI_data> &pointofinterest) {
    for (int featureID = 0; featureID < getNumFeatures(); ++featureID) {
        features[featureID].featurePoints.clear();
    }
    for (int stID = 0; stID < getNumStreets(); ++stID) {
        streets[stID].streetSegments.clear();
    }
    intersections.clear();
    streets.clear();
    features.clear();
    pointofinterest.clear();
}

//remove highlighted intersections
void clearHighlight(std::vector<intersection_data> &intersections, std::vector<POI_data> &pointofinterest) {
    for (int i = 0; i < intersections.size(); ++i) {
        if (intersections[i].highlight) {
            intersections[i].highlight = false;
           
        }
        if (intersections[i].highlight2) {
            intersections[i].highlight2 = false;
           
        }
    }
    for(int i=0;i<pointofinterest.size();++i){
            pointofinterest[i].isAmenity=false;
    }   
}

//returns counts of POI to draw
double determineCount(double worldRatio) {
    if (worldRatio > 50000) return 10;
    if (worldRatio > 20000) return 10;
    if (worldRatio > 7000) return 10;
    if (worldRatio > 400) return 10;
    else return 0;
}

//draws POI
void POIdraw(ezgl::renderer *g, std::vector<POI_data> &pointofinterest, double worldRatio, int &countPOI) {
    ezgl::surface *cafe = g->load_png("libstreetmap/resources/cafe.png");
    ezgl::surface *bar = g->load_png("libstreetmap/resources/bar.png");
    ezgl::surface *pub = g->load_png("libstreetmap/resources/pub.png");
    ezgl::surface *restaurant = g->load_png("libstreetmap/resources/restaurant.png");
    ezgl::surface *fastfood = g->load_png("libstreetmap/resources/fastfood.png");
    ezgl::surface *icecream = g->load_png("libstreetmap/resources/icecream.png");


 for (int numPOI = 0; numPOI < getNumPointsOfInterest(); ++numPOI) {
         float x = x_from_lon(pointofinterest[numPOI].position.lon());
            float y = y_from_lat(pointofinterest[numPOI].position.lat());
            std::string POItype = pointofinterest[numPOI].type;
            ezgl::point2d textPoint = {x, y - 0.000002};
     if(pointofinterest[numPOI].isAmenity){
                g->set_color(ezgl::RED);
                g->fill_arc({x, y}, 0.000005, 0, 360);
            }
 }
    for (int numPOI = 0; numPOI < getNumPointsOfInterest(); ++numPOI) {
        if (g->get_visible_world().contains(p2d_from_ll(pointofinterest[numPOI].position)) && worldRatio > 300 && countPOI < determineCount(worldRatio)) {
            float x = x_from_lon(pointofinterest[numPOI].position.lon());
            float y = y_from_lat(pointofinterest[numPOI].position.lat());
            std::string POItype = pointofinterest[numPOI].type;
            ezgl::point2d textPoint = {x, y - 0.000002};
            
            
            if (pointofinterest[numPOI].isResto && pointofinterest[numPOI].isCafe) {
                g->draw_surface(restaurant,{x, y});
                poiGetName(g, pointofinterest[numPOI].name, POItype, worldRatio, textPoint);
                countPOI++;
            } else if (pointofinterest[numPOI].isResto && pointofinterest[numPOI].isPub) {
                g->draw_surface(restaurant,{x, y});
                poiGetName(g, pointofinterest[numPOI].name, POItype, worldRatio, textPoint);
                countPOI++;
            } else if (pointofinterest[numPOI].isResto) {
                g->draw_surface(restaurant,{x, y});
                poiGetName(g, pointofinterest[numPOI].name, POItype, worldRatio, textPoint);
                countPOI++;
            } else if (pointofinterest[numPOI].isCafe) {
                g->draw_surface(cafe,{x, y});
                poiGetName(g, pointofinterest[numPOI].name, POItype, worldRatio, textPoint);
                countPOI++;
            } else if (pointofinterest[numPOI].isBar) {
                g->draw_surface(bar,{x, y});
                poiGetName(g, pointofinterest[numPOI].name, POItype, worldRatio, textPoint);
                countPOI++;
            } else if (pointofinterest[numPOI].isFastFood) {
                g->draw_surface(fastfood,{x, y});
                poiGetName(g, pointofinterest[numPOI].name, POItype, worldRatio, textPoint);
                countPOI++;
            } else if (pointofinterest[numPOI].isPub) {
                g->draw_surface(pub,{x, y});
                poiGetName(g, pointofinterest[numPOI].name, POItype, worldRatio, textPoint);
                countPOI++;
            } else if (pointofinterest[numPOI].isIce) {
                g->draw_surface(icecream,{x, y});
                poiGetName(g, pointofinterest[numPOI].name, POItype, worldRatio, textPoint);
                countPOI++;

            }


        }

    }
    g->free_surface(cafe);
    g->free_surface(restaurant);
    g->free_surface(pub);
    g->free_surface(bar);
    g->free_surface(fastfood);
    g->free_surface(icecream);

}



//gets the names for POI
void poiGetName(ezgl::renderer *g, std::string poiName, std::string POItype, double worldRatio, ezgl::point2d textPoint) {
    if ((POItype == "restaurant" || POItype == "fast_food" || POItype == "food_court" || POItype == "pub" || POItype == "biergarten" || POItype == "ice_cream" || POItype == "bar" || POItype == "cafe") && worldRatio > 1000) {
        g->set_color(ezgl::BLACK);
        g->format_font("monospace", ezgl::font_slant::normal, ezgl::font_weight::normal, 10);
        g->set_font_size(10);
        g->set_text_rotation(0);
        g->draw_text(textPoint, poiName);
    }
}

//returns specs of features
ezgl::color featureSpecification(int featureType) {
    if (featureType == 1 || featureType == 5 || featureType == 7 || featureType == 8) {
        return ezgl::color(100, 225, 105, 255);
    } else if (featureType == 3 || featureType == 4 || featureType == 9) {
        return ezgl::color(140, 200, 255, 255);
    } else if (featureType == 0 || featureType == 6) {
        return ezgl::color(255, 255, 200, 255);
    } else if (featureType == 2) {
        return ezgl::color(210, 190, 140, 255);
    }
}

//draw arrow for streets
void draw_arrow(ezgl::renderer* g, std::vector<int> segments) {
    for (auto segIdx : segments) {
        if (segIdx % 5 == 0) {
            InfoStreetSegment segInfo = getInfoStreetSegment(segIdx);
            LatLon intFrom = getIntersectionPosition(segInfo.from);
            LatLon intTo = getIntersectionPosition(segInfo.to);
            std::string segType = getSegmentType(segIdx);
            g->set_color(ezgl::GREY_55);
            double angle = atan2(y_from_lat(intTo.lat()) - y_from_lat(intFrom.lat()), (x_from_lon(intTo.lon()) - x_from_lon(intFrom.lon())));
            g->set_text_rotation(angle / DEGREE_TO_RADIAN);
            double width, height;
            if (segType == "motorway" || segType == "primary") {
                width = 0.001;
                height = 0.0001;
            } else if (segType == "secondary" || segType == "tertiary") {
                width = 0.0001;
                height = 0.00001;
            } else {
                width = 0.0005;
                height = 0.00005;
            }
            g->format_font("monospace", ezgl::font_slant::normal, ezgl::font_weight::normal, 12);
            g->draw_text({(x_from_lon(intFrom.lon()) + x_from_lon(intTo.lon())) / 2, (y_from_lat(intFrom.lat()) + y_from_lat(intTo.lat())) / 2}, "->", width, height);
            g->set_text_rotation(0);
        }
    }
}

//call at the end so the streets dont go over it also algorithm was taken online

void featureGetName(ezgl::renderer *g, std::string featureName, FeatureType featureType, double worldRatio, std::vector<ezgl::point2d> points, int fID) {
    std::vector<double> xPoints, yPoints;
    if ((featureType == 1 || featureType == 2 || featureType == 3 || featureType == 4 || featureType == 5) && worldRatio > 400) {
        for (int i = 0; i < points.size(); ++i) {
            xPoints.push_back(points[i].x);
            yPoints.push_back(points[i].y);
        }
        ezgl::point2d p2 = {0, 0};
        double signedA = 0;
        for (int i = 0; i < points.size(); ++i) {
            double x1 = points[i].x;
            double y1 = points[i].y;
            double x2 = points[(i + 1) % points.size()].x;
            double y2 = points[(i + 1) % points.size()].y;
            double a = x1 * y2 - x2*y1;
            signedA += a;
            p2.x += (x1 + x2) * a;
            p2.y += (y1 + y2) * a;

        }
        signedA *= 0.5;
        p2.x /= (6.0 * signedA);
        p2.y /= (6.0 * signedA);

        g->set_color(ezgl::BLACK);
        g->format_font("monospace", ezgl::font_slant::normal, ezgl::font_weight::normal, 10);
        g->set_font_size(9);

        
        if (find_feature_area(fID) < 1000 && worldRatio > 50000) g->draw_text(p2, featureName);
        else if (find_feature_area(fID) < 3000 && worldRatio > 25000) g->draw_text(p2, featureName);
        else if (find_feature_area(fID) < 5000 && worldRatio > 10000) g->draw_text(p2, featureName);
        else if (find_feature_area(fID) < 50000 && worldRatio > 5000) g->draw_text(p2, featureName);
        else if (find_feature_area(fID) > 100000 && worldRatio > 50) g->draw_text(p2, featureName);
    }
}

//average city latitude
double latAvg(double min_lat, double max_lat) {
    lat_avg = (min_lat + max_lat) / 2.0;
    return lat_avg;
}

//some conversion helper functions:

double x_from_lon(double lon) {
    return DEGREE_TO_RADIAN * lon * cos(DEGREE_TO_RADIAN * lat_avg);

}

double y_from_lat(double lat) {
    return DEGREE_TO_RADIAN *lat;
}

double lon_from_x(double x) {
    return x / (DEGREE_TO_RADIAN * cos(DEGREE_TO_RADIAN * lat_avg));
}

double lat_from_y(double y) {
    return y / DEGREE_TO_RADIAN;
}

ezgl::point2d p2d_from_ll(LatLon ll) {
    return
    {
        x_from_lon(ll.lon()), y_from_lat(ll.lat())
    };
}

LatLon ll_from_p2d(ezgl::point2d p2d) {
    return LatLon(lat_from_y(p2d.y), lon_from_x(p2d.x));
}


// Get the centre point of a street segment in x,y coords

ezgl::point2d getCentreOfSeg(LatLon fromPoint, LatLon toPoint) {
    float x1 = x_from_lon(fromPoint.lon());
    float y1 = y_from_lat(fromPoint.lat());
    float x2 = x_from_lon(toPoint.lon());
    float y2 = y_from_lat(toPoint.lat());

    ezgl::point2d seg_centre = {(x1 + x2) / 2, (y1 + y2) / 2};
    return seg_centre;
}
//Get the rotation in degrees of a street segment, in the range [-90, +90] for the purposes
//of setting text rotation so that street names are written left to right and not upside down

float getSegRotation(LatLon fromPoint, LatLon toPoint) {
    float x1 = x_from_lon(fromPoint.lon());
    float y1 = y_from_lat(fromPoint.lat());
    float x2 = x_from_lon(toPoint.lon());
    float y2 = y_from_lat(toPoint.lat());

    float segment_rotation = atan2((y2 - y1), (x2 - x1)) / (M_PI/180.0);
    if (segment_rotation >= -90){
        segment_rotation -= 180;
    }
    if(segment_rotation <= 90){
        segment_rotation += 180;
    }
    return segment_rotation;
}

//get info related to either searched POI or intersection
std::string getInfo(std::vector<intersection_data> &intersections, std::vector<POI_data> &pointofinterest) {
    std::string information;
    for (int i = 0; i < intersections.size(); ++i) {
        if (intersections[i].highlight) {
            information.append("Intersection 1: " + intersections[i].name + "\n\n");
        }
         if (intersections[i].highlight2) {
            information.append("Intersection 2: " + intersections[i].name + "\n\n");
        }
    }
    for (int i = 0; i < pointofinterest.size(); ++i) {
        if (pointofinterest[i].isAmenity) {
            OSMID osmPOI= getPointOfInterestOSMNodeID(i);
            std::map<OSMID, std::string>::iterator it = streetName.find(osmPOI);
             std::map<OSMID, std::string>::iterator it2 = streetNum.find(osmPOI);
             if(((*it2).second).empty()==false){
            information.append(pointofinterest[i].name  + ": " + (*it2).second + " " + (*it).second + "\n");
             }
            
        }
    }
    return information;
}

//for finding directions it returns the clicked on intersection names
std::pair<std::string,std::string> getNameInfo(std::vector<intersection_data> &intersections) {
    std::string information;
    std::pair<std::string,std::string> fromTo;
    for (int i = 0; i < intersections.size(); ++i) {
        if (intersections[i].highlight) {
           fromTo.first= replaceAND(intersections[i].name);
        }
         if (intersections[i].highlight2) {
            fromTo.second=replaceAND(intersections[i].name) ;
        }
    }

    if (fromTo.first.empty())  fromTo.first= "From Intersection  ";
    if (fromTo.second.empty())  fromTo.second= "To Intersection  ";
    return fromTo;
}

//replace & in names with - for searching purposes
std::string replaceAND(std::string intersectionName){
    for(int i=0;i<intersectionName.size();++i){
        if(intersectionName[i]=='&'){
            intersectionName[i]='-';
            return intersectionName;
        }
    }
    return intersectionName;
   
}
