///*
// * To change this license header, choose License Headers in Project Properties.
// * To change this template file, choose Tools | Templates
// * and open the template in the editor.
// */
//
//#include "m2.h"
//#include "m1.h"
//#include "draw.h"
//#include "point.hpp"
//#include "StreetsDatabaseAPI.h"
//#include "OSMDatabaseAPI.h"
//#include "ezgl/application.hpp"
//#include "ezgl/graphics.hpp"
//#include "ezgl/rectangle.hpp"
//#include <cmath>
//
//
//void draw_arrow(ezgl::renderer* g, std::vector<int> segments) {
//    for (auto segIdx : segments) {
//        if (segIdx % 5 == 0) {
//            InfoStreetSegment segInfo = getInfoStreetSegment(segIdx);
//            LatLon intFrom = getIntersectionPosition(segInfo.from);
//            LatLon intTo = getIntersectionPosition(segInfo.to);
//            std::string segType = getSegmentType(segIdx);
//            g->set_color(ezgl::GREY_55);
//            double angle = atan2(y_from_lat(intTo.lat()) - y_from_lat(intFrom.lat()), (x_from_lon(intTo.lon()) - x_from_lon(intFrom.lon())));
//            g->set_text_rotation(angle / DEGREE_TO_RADIAN);
//            double width, height;
//            if (segType == "motorway" || segType == "primary") {
//                width = 0.001;
//                height = 0.0001;
//            } else if (segType == "secondary" || segType == "tertiary") {
//                width = 0.0001;
//                height = 0.00001;
//            } else {
//                width = 0.0005;
//                height = 0.00005;
//            }
//            g->format_font("monospace", ezgl::font_slant::normal, ezgl::font_weight::normal, 12);
//            g->draw_text({(x_from_lon(intFrom.lon()) + x_from_lon(intTo.lon())) / 2, (y_from_lat(intFrom.lat()) + y_from_lat(intTo.lat())) / 2}, "->", width, height);
//            g->set_text_rotation(0);
//        }
//    }
//}
//
//
//void drawStreet(ezgl::renderer *g, std::vector<int> segments) {
//    for (auto segID : segments) {   
//
//        InfoStreetSegment segInfo = getInfoStreetSegment(segID);
//        int segID2=segID;
//        LatLon intFrom = getIntersectionPosition(segInfo.from);
//        LatLon intTo = getIntersectionPosition(segInfo.to);
//      
//        if (g->get_visible_world().contains(p2d_from_ll(intFrom)) || g->get_visible_world().contains(p2d_from_ll(intTo))) {
//            if (segInfo.curvePointCount == 0) {
//                float x1 = x_from_lon(intFrom.lon());
//                float y1 = y_from_lat(intFrom.lat());
//                float x2 = x_from_lon(intTo.lon());
//                float y2 = y_from_lat(intTo.lat());
//
//                g->draw_line({x1, y1},
//                {
//                    x2, y2 });
//
//            } else {
//                for (int i = 0; i < segInfo.curvePointCount - 1; i++) {
//                    LatLon curveP = getStreetSegmentCurvePoint(i, segID);
//                    LatLon curveP2 = getStreetSegmentCurvePoint(i + 1, segID);
//                    float x1 = x_from_lon(curveP.lon());
//                    float y1 = y_from_lat(curveP.lat());
//                    float x2 = x_from_lon(curveP2.lon());
//                    float y2 = y_from_lat(curveP2.lat());
//                    g->draw_line({x1, y1},
//                    {
//                        x2, y2 });
//                    //                if(segID % 5 ==0){
//                    //                    draw_arrow(g, segInfo.oneWay, intFrom, intTo);
//                    //                }
//
//                }
//                LatLon curveFirst = getStreetSegmentCurvePoint(0, segID);
//                LatLon curveLast = getStreetSegmentCurvePoint(segInfo.curvePointCount - 1, segID);
//                float x1 = x_from_lon(curveFirst.lon());
//                float y1 = y_from_lat(curveFirst.lat());
//                float x2 = x_from_lon(intFrom.lon());
//                float y2 = y_from_lat(intFrom.lat());
//                g->draw_line({x1, y1},
//                {
//                    x2, y2 });
//
//                float x3 = x_from_lon(curveLast.lon());
//                float y3 = y_from_lat(curveLast.lat());
//                float x4 = x_from_lon(intTo.lon());
//                float y4 = y_from_lat(intTo.lat());
//                g->draw_line({x3, y3},
//                {
//                    x4, y4 });
//
//            }
//        }
//      
//    }
//}
//
//std::pair<int,ezgl::color> segmentSpecification(std::string segmentType) {
//    if (segmentType == "motorway") {
//        return std::make_pair(6, ezgl::color(255, 140, 100, 255));
//    } else if (segmentType == "trunk") {
//        return std::make_pair(6, ezgl::PLUM);
//    } else if (segmentType == "primary") {
//        return std::make_pair(5, ezgl::ORANGE);
//    } else if (segmentType == "secondary") {
//        return std::make_pair(6, ezgl::BISQUE);
//    } else if (segmentType == "tertiary") {
//        return std::make_pair(5, ezgl::BISQUE);
//    } else if (segmentType == "unclassified") {
//        return std::make_pair(5, ezgl::WHITE);
//    } else if (segmentType == "residential") {
//        return std::make_pair(5, ezgl::WHITE);
//    }  
//}
//
//void setColourandWidth(ezgl::renderer *g, std::string segType) {
//
//    std::pair<double, ezgl::color> WidthandColor = segmentSpecification(segType);
//    g->set_color(WidthandColor.second);
//    g->set_line_width(WidthandColor.first);
//
//}
//
