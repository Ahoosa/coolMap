
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "m3.h"
#include "m2.h"
#include "m1.h"
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
#include <stdlib.h>
#include "LatLon.h"



void draw_map_blank_canvas();
int highlightCount = 0;
void FastS(GtkWidget *widget, ezgl::application *app);
void BarS(GtkWidget *widget, ezgl::application *app);
void CafeS(GtkWidget *widget, ezgl::application *app);
void RestoS(GtkWidget *widget, ezgl::application *app);
void IceS(GtkWidget *widget, ezgl::application *app);
void PubS(GtkWidget *widget, ezgl::application *app);
void AnyArrow(GtkWidget *widget, ezgl::application *app);
void walk(GtkWidget *widget, ezgl::application *app);
int countPOI;
void zoomTofit(ezgl::application *application, LatLon zoompoint);
void ZoomButton(GtkWidget *widget, ezgl::application *app);
double worldRatio;
void initial_setup(ezgl::application *application, bool new_window);
void act_on_mouse_button(ezgl::application *app, GdkEventButton *event, double x, double y);
void act_on_mouse_click(ezgl::application* app, GdkEventButton* event, double x, double y);
void find_button(GtkWidget *widget, ezgl::application *application);
void draw_main_canvas(ezgl::renderer *g);
double max_lat, min_lat, max_lon, min_lon;
void zoomIn(ezgl::application *application, LatLon zoompoint);
void zoomOut(ezgl::application *application);
void city_selecter(GtkWidget *widget, ezgl::application *application);
void pop_up(GtkWidget *widget, ezgl::application *app);
void help_pop_up(GtkWidget *widget, ezgl::application *app);
void find_directions(GtkWidget *widget, ezgl::application *app, ezgl::renderer* g);
void drawPath(ezgl::renderer *g, std::vector<IntersectionIndex> path, IntersectionIndex fromID, IntersectionIndex toID, bool walk);
void giveDirections(std::vector<StreetSegmentIndex> path);
void giveDirectionsWalking(std::vector<IntersectionIndex> walk, std::vector<IntersectionIndex> drive);
void drawPathWalking(ezgl::renderer *g, std::vector<IntersectionIndex> path);

bool startDisplay;
bool searchEnabled;
bool walkEnabled = false;
double walkTime, walkSpeed;
void act_on_mouse_move(ezgl::application* app,
        GdkEventButton* event,
        double x, double y);
void on_dialog_response(GtkDialog *dialog, gint response_id, gpointer user_data);
void search_dialog_response(GtkDialog *dialog, gint response_id, gpointer user_data);
void error_dialog_response(GtkDialog *dialog, gint response_id, gpointer user_data);
void draw_arrow(ezgl::renderer* g, bool one_way, LatLon from, LatLon to);
void featureNames(ezgl::renderer* g);
bool sortbysec(const std::pair<int, double> &a, const std::pair<int, double> &b);

GtkEntry* entry1;
GtkEntry* entry2;
GtkEntry* entry6;

std::vector<int> segMotor, segTrunk, segPrim, segSec, segTer, segUnc, segRes;
std::vector<intersection_data> intersections;
std::vector<feature_data> features;
std::vector<street_data> streets;
std::vector<POI_data> pointofinterest;
std::pair<IntersectionIndex, IntersectionIndex> directionIntersections;

void draw_map() {

    std::pair<std::pair<double, double>, std::pair<double, double>> MaxandMin = getMaxMin();

    max_lat = MaxandMin.first.first;
    max_lon = MaxandMin.first.second;
    min_lat = MaxandMin.second.first;
    min_lon = MaxandMin.second.second;
    double lat_avg = latAvg(min_lat, max_lat);


    fillVec(intersections, features, streets, pointofinterest);

    ezgl::application::settings settings;
    settings.main_ui_resource = "libstreetmap/resources/main.ui";
    settings.window_identifier = "MainWindow";
    settings.canvas_identifier = "MainCanvas";
    ezgl::application application(settings);


    ezgl::rectangle initial_world({x_from_lon(min_lon), y_from_lat(min_lat)},
    {
        x_from_lon(max_lon), y_from_lat(max_lat) });
    application.add_canvas("MainCanvas", draw_main_canvas, initial_world);

    ezgl::rectangle changed_world = initial_world;
    worldRatio = initial_world.area() / changed_world.area();
    startDisplay = false;


    application.run(initial_setup, act_on_mouse_click, act_on_mouse_move, nullptr);

}

void initial_setup(ezgl::application *application, bool /*new_window*/) {

    GObject *FindButton = application->get_object("FindButton");
    g_signal_connect(G_OBJECT(FindButton), "clicked", G_CALLBACK(find_button), application);

    GtkComboBox *changedCity = (GtkComboBox *) application->get_object("CityList");
    g_signal_connect(G_OBJECT(changedCity), "changed", G_CALLBACK(city_selecter), application);

    GObject *DisplayB = application->get_object("Display");
    g_signal_connect(G_OBJECT(DisplayB), "clicked", G_CALLBACK(pop_up), application);

    GObject *ZoomInButton = application->get_object("ZoomInButton");
    g_signal_connect(G_OBJECT(ZoomInButton), "clicked", G_CALLBACK(ZoomButton), application);

    GObject *ZoomOutButton = application->get_object("ZoomOutButton");
    g_signal_connect(G_OBJECT(ZoomOutButton), "clicked", G_CALLBACK(ZoomButton), application);

    GObject *ZoomFitButton = application->get_object("ZoomFitButton");
    g_signal_connect(G_OBJECT(ZoomFitButton), "clicked", G_CALLBACK(ZoomButton), application);

    GObject *DirectionsButton = application->get_object("DirectionsButton");
    g_signal_connect(G_OBJECT(DirectionsButton), "clicked", G_CALLBACK(find_directions), application);

    GObject *HelpButton = application->get_object("HelpButton");
    g_signal_connect(G_OBJECT(HelpButton), "clicked", G_CALLBACK(help_pop_up), application);

    GObject *WalkDrive = application->get_object("WalkDrive");
    g_signal_connect(G_OBJECT(WalkDrive), "clicked", G_CALLBACK(walk), application);

    GObject *RestoB = application->get_object("RestoB");
    g_signal_connect(G_OBJECT(RestoB), "toggled", G_CALLBACK(RestoS), application);
    GObject *BarB = application->get_object("BarB");
    g_signal_connect(G_OBJECT(BarB), "toggled", G_CALLBACK(BarS), application);
    GObject *CafeB = application->get_object("CafeB");
    g_signal_connect(G_OBJECT(CafeB), "toggled", G_CALLBACK(CafeS), application);
    GObject *PubB = application->get_object("PubB");
    g_signal_connect(G_OBJECT(PubB), "toggled", G_CALLBACK(PubS), application);
    GObject *IceB = application->get_object("IceB");
    g_signal_connect(G_OBJECT(IceB), "toggled", G_CALLBACK(IceS), application);
    GObject *FastB = application->get_object("FastB");
    g_signal_connect(G_OBJECT(FastB), "toggled", G_CALLBACK(FastS), application);

}

void draw_main_canvas(ezgl::renderer *g) {

    g->draw_rectangle({x_from_lon(min_lon), y_from_lat(min_lat)},
    {
        x_from_lon(max_lon), y_from_lat(max_lat)
    });
    g->set_color(ezgl::color(236, 230, 222, 200));
    g->fill_rectangle(g->get_visible_world());




    std::vector<std::pair<int, double>> indexArea;
    indexArea.resize(getNumFeatures());

    for (int featureIndex = 0; featureIndex < getNumFeatures(); ++featureIndex) {
        indexArea[featureIndex].first = featureIndex;
        indexArea[featureIndex].second = find_feature_area(featureIndex);
    }
    std::sort(indexArea.begin(), indexArea.end(), sortbysec);

    for (int j = 0; j < features.size(); ++j) {
        int index = indexArea[j].first;
        if (features[j].closed_feature && features[j].featurePoints.size() > 1) {
            g->set_color(featureSpecification(features[j].type));

            g->fill_poly(features[j].featurePoints);
        } else {
            //go through every feature, then go through every pt of that feature & get the x and y and draw it
            for (int point = 0; point < getFeaturePointCount(index); ++point) {
                float x1 = (features[j].featurePoints[point].x);
                float y1 = (features[j].featurePoints[point].y);

                float width1 = 0.00001;

                if (g->get_visible_world().contains(features[j].featurePoints[point])) {
                    g->set_color(featureSpecification(features[j].type));
                    g->fill_arc({x1, y1}, 0.0000005, 0, 360);
                }

            }
        }
    }





    for (int stID = 0; stID < streets.size(); ++stID) {


        std::string segmentType = streets[stID].type;
        if (segmentType == "motorway") segMotor.push_back(stID);
        if (segmentType == "trunk") segTrunk.push_back(stID);
        if (segmentType == "primary") segPrim.push_back(stID);
        if (segmentType == "secondary") segSec.push_back(stID);
        if (segmentType == "tertiary") segTer.push_back(stID);
        if (segmentType == "unclassified") segUnc.push_back(stID);
        if (segmentType == "residential") segRes.push_back(stID);

    }



    if (worldRatio > 400) {

        for (auto stID : segRes) {
            std::string segmentType = "residential";
            setColourandWidth(g, segmentType);
            drawStreet(g, streets[stID].streetSegments);
            if (worldRatio > 5000) draw_arrow(g, streets[stID].streetSegments);
        }
        for (auto stID : segUnc) {
            std::string segmentType = "unclassified";
            setColourandWidth(g, segmentType);
            drawStreet(g, streets[stID].streetSegments);
            if (worldRatio > 5000) draw_arrow(g, streets[stID].streetSegments);
        }

    }
    if (worldRatio > 50) {
        for (auto stID : segTer) {
            std::string segmentType = "tertiary";
            setColourandWidth(g, segmentType);
            drawStreet(g, streets[stID].streetSegments);
            if (worldRatio > 5000) draw_arrow(g, streets[stID].streetSegments);
        }

    }
    if (worldRatio > 20) {
        for (auto stID : segSec) {
            std::string segmentType = "secondary";
            setColourandWidth(g, segmentType);
            drawStreet(g, streets[stID].streetSegments);
            if (worldRatio > 5000) draw_arrow(g, streets[stID].streetSegments);
        }

        for (auto stID : segPrim) {
            std::string segmentType = "primary";
            setColourandWidth(g, segmentType);
            drawStreet(g, streets[stID].streetSegments);
            if (worldRatio > 5000) draw_arrow(g, streets[stID].streetSegments);
        }


    }

    for (auto stID : segTrunk) {
        std::string segmentType = "trunk";
        setColourandWidth(g, segmentType);
        drawStreet(g, streets[stID].streetSegments);
        //      if (worldRatio > 3000)  draw_arrow(g, streets[stID].streetSegments);
    }
    for (auto stID : segMotor) {
        std::string segmentType = "motorway";
        setColourandWidth(g, segmentType);
        drawStreet(g, streets[stID].streetSegments);
        //       if (worldRatio > 3000) draw_arrow(g, streets[stID].streetSegments);
    }





    segMotor.clear();
    segTrunk.clear();
    segPrim.clear();
    segSec.clear();
    segTer.clear();
    segUnc.clear();
    segRes.clear();

    for (int j = 0; j < features.size(); ++j) {
        if (features[j].name != "<noname>") {
            featureGetName(g, features[j].name, features[j].type, worldRatio, features[j].featurePoints, j);
        }
    }



    countPOI = 0;
    POIdraw(g, pointofinterest, worldRatio, countPOI);


    for (int i = 0; i < intersections.size(); ++i) {
        if (g->get_visible_world().contains(p2d_from_ll(intersections[i].position))) {
            float x = x_from_lon(intersections[i].position.lon());
            float y = y_from_lat(intersections[i].position.lat());


            if (intersections[i].highlight) {
                g->set_color(ezgl::RED);
                g->fill_arc({x, y}, 0.000005, 0, 360);
            }
            if (intersections[i].highlight2) {
                g->set_color(ezgl::RED);
                g->fill_arc({x, y}, 0.000005, 0, 360);
            }
        }
    }
    if (searchEnabled) {
        if (walkEnabled) {
            std::vector<StreetSegmentIndex> walk = find_path_with_walk_to_pick_up(directionIntersections.first, directionIntersections.second, 1.0, walkSpeed, walkTime).first;
            std::vector<StreetSegmentIndex> drive = find_path_with_walk_to_pick_up(directionIntersections.first, directionIntersections.second, 1.0, walkSpeed, walkTime).second;
            if (walk.size() != 0) {
                drawPathWalking(g, walk);
            }
            if (drive.size() != 0) {
                drawPath(g, drive, getInfoStreetSegment(drive.front()).from, getInfoStreetSegment(drive.back()).to, walkEnabled);
            }
        }
        else {
            drawPath(g, find_path_between_intersections(directionIntersections.first, directionIntersections.second, 1.0), directionIntersections.first, directionIntersections.second, walkEnabled);
        }
    }
}

void draw_map_blank_canvas() {
    ezgl::application::settings settings;
    settings.main_ui_resource = "libstreetmap/resourcea/main.ui";
    settings.window_identifier = "MainWindow";
    settings.canvas_identifier = "MainCanvas";

    ezgl::application application(settings);

    ezgl::rectangle initial_world({x_from_lon(min_lon), y_from_lat(min_lat)},
    {
        x_from_lon(max_lon), y_from_lat(max_lat)
    });
    application.add_canvas("MainCanvas", draw_main_canvas, initial_world);

    application.run(nullptr, nullptr, nullptr, nullptr);
}

void act_on_mouse_click(ezgl::application* app,
        GdkEventButton* event,
        double x, double y) {
    zoomOut(app);

    std::cout << "Mouse clicked at (" << x << "," << y << ")\n";
    LatLon position = LatLon(lat_from_y(y), lon_from_x(x));
    int id = find_closest_intersection(position);
    std::cout << "Closest Intersection: "
            << intersections[id].name << "\n";

    if (highlightCount == 2) {
        clearHighlight(intersections, pointofinterest);
        highlightCount = 0;
    }
    if (intersections[id].highlight == false && highlightCount == 0) {
        intersections[id].highlight = true;
        highlightCount++;
    } else if (highlightCount == 1) {
        intersections[id].highlight2 = true;
        highlightCount++;
    }

    zoomIn(app, intersections[id].position);
    worldRatio = app->get_canvas(app->get_main_canvas_id())->get_camera().get_initial_world().area() / app->get_canvas(app->get_main_canvas_id())->get_camera().get_world().area();
    //force a refresh (redraw) after mouse click
    app->refresh_drawing();
}

void act_on_mouse_move(ezgl::application* app,
        GdkEventButton* event,
        double x, double y) {
    if (startDisplay == false && worldRatio > 3000) {

        ezgl::point2d hover = {x, y};

        for (int i = 0; i < streets.size(); ++i) {
            showSegNames(app, streets[i].streetSegments, hover, startDisplay);
        }
    }

}

void find_button(GtkWidget *widget, ezgl::application *application) {

    GtkEntry* text_entry = (GtkEntry *) application->get_object("TextInput");
    zoomOut(application);
    std::string text = gtk_entry_get_text(text_entry);
    std::stringstream iss(text);
    std::string a;
    bool foundIt = false;
    std::vector<int> stID1;
    std::vector<int> stID2;
    clearHighlight(intersections, pointofinterest);
    int c = 0;

    for (int i = 0; i < pointofinterest.size(); ++i) {
        boost::algorithm::to_upper(pointofinterest[i].name);
        boost::algorithm::to_upper(text);
        if (text == pointofinterest[i].name) {
            std::cout << pointofinterest[i].name << std::endl;
            pointofinterest[i].isAmenity = true;
            foundIt = true;
            if (c == 0) {
               zoomOut(application);
                zoomIn(application, pointofinterest[i].position);
                c++;
            }
        }
    }
    if (foundIt == false) {
        while (getline(iss, a, '-')) {

            if (stID1.size() == 0) {
                stID1 = find_street_ids_from_partial_street_name(a);
            }
            if (stID1.size() != 0) {
                stID2 = find_street_ids_from_partial_street_name(a);
            }
        }


        std::vector<int> intOfTwo;
        int m = 0;
        for (auto i : stID1) {
            for (auto j : stID2) {
                intOfTwo = find_intersections_of_two_streets(std::make_pair(i, j));
                if (intOfTwo.size() != 0) {
                    for (auto id : intOfTwo) {

                        if (m == 0) {
                            zoomOut(application);
                            zoomIn(application, intersections[id].position);
                            m++;
                        }
                        intersections[id].highlight = true;
                    }
                }
                intOfTwo.clear();
            }
        }

    }
    worldRatio = application->get_canvas(application->get_main_canvas_id())->get_camera().get_initial_world().area() / application->get_canvas(application->get_main_canvas_id())->get_camera().get_world().area();
    application->refresh_drawing();

}

void zoomIn(ezgl::application *application, LatLon zoompoint) {
    ezgl::zoom_in(application->get_canvas("MainCanvas"), application->get_canvas("MainCanvas")->get_camera().world_to_screen(p2d_from_ll(zoompoint)), 18);
}

void zoomOut(ezgl::application *application) {
    ezgl::zoom_fit(application->get_canvas("MainCanvas"), application->get_canvas("MainCanvas")->get_camera().get_initial_world());
    worldRatio = application->get_canvas(application->get_main_canvas_id())->get_camera().get_initial_world().area() / application->get_canvas(application->get_main_canvas_id())->get_camera().get_world().area();
}

void zoomTofit(ezgl::application *application, LatLon zoompoint) {
    ezgl::zoom_in(application->get_canvas("MainCanvas"), application->get_canvas("MainCanvas")->get_camera().world_to_screen(p2d_from_ll(zoompoint)), 2);
}

void city_selecter(GtkWidget *widget, ezgl::application *app) {
    GtkComboBoxText* text_entry = (GtkComboBoxText *) app->get_object("CityList");
    std::string cityName = gtk_combo_box_text_get_active_text(text_entry);
    std::cout << cityName << std::endl;


    clearData(intersections, features, streets, pointofinterest);
    close_map();

    bool load_success = load_map(find_city(cityName));
    if (load_success) {
        std::cout << "load successful" << std::endl;
    }

    fillVec(intersections, features, streets, pointofinterest);
    app->update_message("Map: " + cityName);
    auto MaxandMin = getMaxMin();
    max_lat = MaxandMin.first.first;
    max_lon = MaxandMin.first.second;
    min_lat = MaxandMin.second.first;
    min_lon = MaxandMin.second.second;

    const ezgl::rectangle new_world(p2d_from_ll(LatLon(max_lat, max_lon)), p2d_from_ll(LatLon(min_lat, min_lon)));
    app->change_canvas_world_coordinates(app->get_main_canvas_id(), new_world);

    worldRatio = 1;
    if (app->get_canvas("MainCanvas")->height() > new_world.height() || app->get_canvas("MainCanvas")->height() > new_world.height()) {
        zoomTofit(app, ll_from_p2d(new_world.center()));
    }

    app->refresh_drawing();

}

void on_dialog_response(GtkDialog *dialog, gint response_id, gpointer user_data) {
    switch (response_id) {
        case GTK_RESPONSE_ACCEPT:
            //get 
            break;
        case GTK_RESPONSE_DELETE_EVENT:
            break;
        case GTK_RESPONSE_REJECT:
            break;
        default:
            break;
    }
    gtk_widget_destroy(GTK_WIDGET(dialog));
}

void search_dialog_response(GtkDialog *dialog, gint response_id, gpointer user_data) {
    ezgl::application *app = (ezgl::application*)user_data;

    switch (response_id) {
        case GTK_RESPONSE_ACCEPT:
        {
            std::string searchIntersection1 = gtk_entry_get_text(entry1);
            std::string searchIntersection2 = gtk_entry_get_text(entry2);

            std::stringstream get(searchIntersection1);
            std::stringstream get2(searchIntersection2);

            std::vector<int> stID1;
            std::vector<int> stID2;
            std::vector<int> stID3;
            std::vector<int> stID4;

            std::string a, b;
            while (getline(get, a, '-')) {

                if (stID1.size() == 0) {
                    stID1 = find_street_ids_from_partial_street_name(a);
                }
                if (stID1.size() != 0) {
                    stID2 = find_street_ids_from_partial_street_name(a);
                }
            }
            while (getline(get2, b, '-')) {

                if (stID3.size() == 0) {
                    stID3 = find_street_ids_from_partial_street_name(b);
                }
                if (stID3.size() != 0) {
                    stID4 = find_street_ids_from_partial_street_name(b);
                }
            }
            if (walkEnabled) {
                std::string walkTimeSpeed = gtk_entry_get_text(entry6);
                std::stringstream get3(walkTimeSpeed);
                std::vector<std::string> Time, Speed;
                std::string d;
                bool loop = false;
                while (getline(get3, d, '-')) {
                    if (!loop) {
                        walkTime = atof(d.c_str())*60;
                        loop = true;
                    }
                    if (loop) {
                        walkSpeed = atof(d.c_str());
                    }
                }
            }

            if (stID1.size() == 0 || stID2.size() == 0 || stID3.size() == 0 || stID4.size() == 0) {
                searchEnabled = false;
                GObject *window;
                GtkWidget *content_area;
                GtkWidget *label;
                GtkWidget* dialog_new;
                GtkWidget* dialog;

                window = app->get_object(app->get_main_window_id().c_str());
                dialog_new = gtk_dialog_new_with_buttons(("Error:"),
                        (GtkWindow*) window,
                        GTK_DIALOG_MODAL,
                        ("Close"),
                        GTK_RESPONSE_REJECT,
                        NULL);
                content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog_new));
                label = gtk_label_new("One of the street names you entered was invalid");
                gtk_container_add(GTK_CONTAINER(content_area), label);
                gtk_widget_show_all(dialog_new);
                g_signal_connect(GTK_DIALOG(dialog_new), "response", G_CALLBACK(error_dialog_response), app);
            } else {
                std::vector<int> intOfTwo1, intOfTwo2;
                int c = 0;
                for (auto i : stID1) {
                    for (auto j : stID2) {
                        intOfTwo1 = find_intersections_of_two_streets(std::make_pair(i, j));
                        if (intOfTwo1.size() != 0) {
                            for (auto iD : intOfTwo1) {
                                directionIntersections.first = iD;
                                if(c==0){
                                   zoomIn(app, intersections[iD].position); 
                                   c++;
                                }      
                                if (c == 0) {
                                    zoomOut(app);
                                    directionIntersections.first = iD;
                                    zoomIn(app, getIntersectionPosition(iD));
                                    c++;
                                }
                            }
                        }
                        intOfTwo1.clear();
                    }
                }

                for (auto k : stID3) {
                    for (auto l : stID4) {
                        intOfTwo2 = find_intersections_of_two_streets(std::make_pair(k, l));
                        if (intOfTwo2.size() != 0) {
                            for (auto id : intOfTwo2) {

                                directionIntersections.second = id;
                            }
                        }
                        intOfTwo2.clear();
                    }
                }

                if (walkEnabled) {
                    std::vector<IntersectionIndex> walk = find_path_with_walk_to_pick_up(directionIntersections.first, directionIntersections.second, 0.0, walkSpeed , walkTime).first;
                    std::vector<IntersectionIndex> drive = find_path_with_walk_to_pick_up(directionIntersections.first, directionIntersections.second, 0.0, walkSpeed , walkTime).second;
                    if (walk.size() == 0 && drive.size() == 0) {
                        searchEnabled = false;

                        GObject *window = app->get_object(app->get_main_window_id().c_str());
                        GtkWidget *dialog1 = gtk_dialog_new_with_buttons(("Error:"),
                                (GtkWindow*) window,
                                GTK_DIALOG_MODAL,
                                ("Close"),
                                GTK_RESPONSE_REJECT,
                                NULL);
                        GtkWidget *content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog1));
                        GtkWidget *label = gtk_label_new("A path could not be found between these intersections");
                        gtk_container_add(GTK_CONTAINER(content_area), label);

                        gtk_widget_show_all(dialog1);
                        g_signal_connect(GTK_DIALOG(dialog), "response", G_CALLBACK(error_dialog_response), app);
                        app->refresh_drawing();
                    } else {
                        searchEnabled = true;
                        giveDirectionsWalking(walk, drive);   
                    } 
                }
                else {
                    std::vector<int> path = find_path_between_intersections(directionIntersections.first, directionIntersections.second, 0.0);
                    if (path.size() == 0) {
                        searchEnabled = false;
                        GObject *window = app->get_object(app->get_main_window_id().c_str());
                        GtkWidget *dialog2 = gtk_dialog_new_with_buttons(("Error:"),
                                (GtkWindow*) window,
                                GTK_DIALOG_MODAL,
                                ("Close"),
                                GTK_RESPONSE_REJECT,
                                NULL);
                        GtkWidget *content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog2));
                        GtkWidget *label = gtk_label_new("A path could not be found between these intersections");
                        gtk_container_add(GTK_CONTAINER(content_area), label);
                        gtk_widget_show_all(dialog2);
                        g_signal_connect(GTK_DIALOG(dialog2), "response", G_CALLBACK(error_dialog_response), app);
                    } else {
                        giveDirections(path);
                        searchEnabled = true;
                    }         
                        app->refresh_drawing();
                }
            }
            break;
        }
        case GTK_RESPONSE_DELETE_EVENT:
            searchEnabled = false;
            break;
        case GTK_RESPONSE_REJECT:
            searchEnabled = false;
            break;
        default:
            break;
    }
    gtk_widget_destroy(GTK_WIDGET(dialog));
}

void error_dialog_response(GtkDialog *dialog, gint response_id, gpointer user_data) {
    std::cout << "response is ";
    switch (response_id) {
            std::cout << "Error in Search\n";
        case GTK_RESPONSE_ACCEPT:
        {
            //get 
            break;
        }
        case GTK_RESPONSE_DELETE_EVENT:
            break;
        case GTK_RESPONSE_REJECT:
            break;
        default:
            break;
    }
    gtk_widget_destroy(GTK_WIDGET(dialog));
}

void find_directions(GtkWidget *widget, ezgl::application *app, ezgl::renderer* g) {
    walkEnabled = false;
    GObject *window;
    GtkWidget *content_area;
    GtkWidget *label;
    GtkWidget* dialog;
    GtkWidget* entry3;
    GtkWidget* entry4;

    window = app->get_object(app->get_main_window_id().c_str());
    dialog = gtk_dialog_new_with_buttons(("Enter Intersections"),
            (GtkWindow*) window,
            GTK_DIALOG_MODAL,
            ("Search"),
            GTK_RESPONSE_ACCEPT,
            ("Close"),
            GTK_RESPONSE_REJECT,
            NULL);
    content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
    label = gtk_label_new("Enter the names of the intersections you wish to travel between");
    gtk_container_add(GTK_CONTAINER(content_area), label);
    entry3 = gtk_entry_new();
    entry4 = gtk_entry_new();
    if (highlightCount > 0) {
        std::pair<std::string, std::string> fromTo = getNameInfo(intersections);
        gtk_entry_set_text((GtkEntry*) entry3, fromTo.first.c_str());
        if (highlightCount == 2) gtk_entry_set_text((GtkEntry*) entry4, fromTo.second.c_str());
    } else {
        gtk_entry_set_text((GtkEntry*) entry3, "From Intersection");
        gtk_entry_set_text((GtkEntry*) entry4, "To Intersection");
    }
    gtk_container_add(GTK_CONTAINER(content_area), entry3);
    gtk_container_add(GTK_CONTAINER(content_area), entry4);
    gtk_widget_show_all(dialog);

    entry1 = (GtkEntry*) entry3;
    entry2 = (GtkEntry*) entry4;

    g_signal_connect(GTK_DIALOG(dialog), "response", G_CALLBACK(search_dialog_response), app);
    app->refresh_drawing();
}

void help_pop_up(GtkWidget *widget, ezgl::application *app) {
    GObject *window;
    GtkWidget *content_area;
    GtkWidget *label;
    GtkWidget *dialog;

    window = app->get_object(app->get_main_window_id().c_str());
    dialog = gtk_dialog_new_with_buttons(("How to use our Cool Map"),
            (GtkWindow*) window,
            GTK_DIALOG_MODAL,
            ("OK"),
            GTK_RESPONSE_ACCEPT,
            NULL);
    content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
    label = gtk_label_new("ARROWS: move around on map interface\nZOOM IN: zoom in and see more details and points of interest.\n\n"
            "ZOOM OUT: zoom out to see fewer details and larger area.\n\n"
            "ZOOM FIT: Fit map to see full rendering\n\nPROCEED: Close map.\n\n"
            "DISPLAY: Give information about intersections and points of interest that you have clicked on or searched for.\n\n"
            "FIND DIRECTIONS: Enter two intersections in format street1-street2, street3-street4 to find directions between them\nDirections will be given for driving only.\n\n"
            "AMENITIES LAYERS: Select which type of amenities you would like to see on your map. Note that amenity pop-ups only appear once you have zoomed in a certain amount.\n\n"
            "SEARCH: search for intersections, streets and amenities by name.\n\n"
            "CLICK AND ZOOM: Click anywhere on the map to find the closest intersection. Map will zoom into the intersection.\n\n"
            "SELECT MAP: Use drop down menu at the bottom to select which map you will load.\n\n"
            "FIND PICKUP: Allows you to use the 'Uber pool' function to find an intersection you can walk to. Must specify walking speed and maximum walking time.\n"
            "The walking path will displayed in pink, and the driving path will be displayed in blue.\n"
            "Intersections must be given in the format 'street1-street2', 'street3-street4', and walking time and speed must given in the form 'max walking time-walking speed'\n\n"
            "CLICK FOR DIRECTIONS: Select any two points on the map to be given the directions between the closest intersections to these points. Can be done with either driving only or driving and walking.\n"
            "Simply click on your two points, which will be highlighted in red, then select either the Find Directions button or the Find Pickup button\n\n"
            );
    gtk_container_add(GTK_CONTAINER(content_area), label);
    gtk_widget_show_all(dialog);
    g_signal_connect(GTK_DIALOG(dialog), "response", G_CALLBACK(on_dialog_response), NULL);
    app->refresh_drawing();
}

void pop_up(GtkWidget *widget, ezgl::application *app) {
    GObject *window;
    GtkWidget *content_area;
    GtkWidget *label;
    GtkWidget *dialog;
    std::string interInformation = getInfo(intersections, pointofinterest);

    window = app->get_object(app->get_main_window_id().c_str());
    dialog = gtk_dialog_new_with_buttons(("Information"),
            (GtkWindow*) window,
            GTK_DIALOG_MODAL,
            ("OK"),
            GTK_RESPONSE_ACCEPT,
            ("Close"),
            GTK_RESPONSE_REJECT,
            NULL);

    content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
    label = gtk_label_new(interInformation.c_str());
    gtk_container_add(GTK_CONTAINER(content_area), label);
    gtk_widget_show_all(dialog);
    g_signal_connect(GTK_DIALOG(dialog), "response", G_CALLBACK(on_dialog_response), NULL);

    app->refresh_drawing();
}

void ZoomButton(GtkWidget *widget, ezgl::application *app) {
    worldRatio = app->get_canvas(app->get_main_canvas_id())->get_camera().get_initial_world().area() / app->get_canvas(app->get_main_canvas_id())->get_camera().get_world().area();
    app->refresh_drawing();
}

void AnyArrow(GtkWidget *widget, ezgl::application *app) {

    app->refresh_drawing();
}

void RestoS(GtkWidget *widget, ezgl::application *app) {
    for (int numPOI = 0; numPOI < getNumPointsOfInterest(); ++numPOI) {
        float x = x_from_lon(pointofinterest[numPOI].position.lon());
        float y = y_from_lat(pointofinterest[numPOI].position.lat());
        std::string POItype = pointofinterest[numPOI].type;


        if (pointofinterest[numPOI].isResto == true) {
            pointofinterest[numPOI].isResto = false;
        } else {
            if (POItype == "restaurant" || POItype == "food_court") {
                pointofinterest[numPOI].isResto = true;
            } else {
                pointofinterest[numPOI].isResto = false;
            }
        }
    }
    app->refresh_drawing();
}

void CafeS(GtkWidget *widget, ezgl::application *app) {
    for (int numPOI = 0; numPOI < getNumPointsOfInterest(); ++numPOI) {

        std::string POItype = pointofinterest[numPOI].type;

        if (pointofinterest[numPOI].isCafe == true) {
            pointofinterest[numPOI].isCafe = false;
        } else {
            if (POItype == "cafe") {
                pointofinterest[numPOI].isCafe = true;
            } else {
                pointofinterest[numPOI].isCafe = false;
            }
        }
    }
    app->refresh_drawing();
}

void BarS(GtkWidget *widget, ezgl::application *app) {
    for (int numPOI = 0; numPOI < getNumPointsOfInterest(); ++numPOI) {
        float x = x_from_lon(pointofinterest[numPOI].position.lon());
        float y = y_from_lat(pointofinterest[numPOI].position.lat());
        std::string POItype = pointofinterest[numPOI].type;


        if (pointofinterest[numPOI].isBar == true) {
            pointofinterest[numPOI].isBar = false;
        } else {
            if (POItype == "bar") {
                pointofinterest[numPOI].isBar = true;
            } else {
                pointofinterest[numPOI].isBar = false;
            }
        }
    }
    app->refresh_drawing();
}

void PubS(GtkWidget *widget, ezgl::application *app) {
    for (int numPOI = 0; numPOI < getNumPointsOfInterest(); ++numPOI) {
        float x = x_from_lon(pointofinterest[numPOI].position.lon());
        float y = y_from_lat(pointofinterest[numPOI].position.lat());
        std::string POItype = pointofinterest[numPOI].type;


        if (pointofinterest[numPOI].isPub == true) {
            pointofinterest[numPOI].isPub = false;
        } else {
            if (POItype == "pub" || POItype == "biergarten") {
                pointofinterest[numPOI].isPub = true;
            } else {
                pointofinterest[numPOI].isPub = false;
            }
        }
    }
    app->refresh_drawing();
}

void IceS(GtkWidget *widget, ezgl::application *app) {
    for (int numPOI = 0; numPOI < getNumPointsOfInterest(); ++numPOI) {
        float x = x_from_lon(pointofinterest[numPOI].position.lon());
        float y = y_from_lat(pointofinterest[numPOI].position.lat());
        std::string POItype = pointofinterest[numPOI].type;


        if (pointofinterest[numPOI].isIce == true) {
            pointofinterest[numPOI].isIce = false;
        } else {
            if (POItype == "ice_cream") {
                pointofinterest[numPOI].isIce = true;
            } else {
                pointofinterest[numPOI].isIce = false;
            }
        }
    }
    app->refresh_drawing();
}

void FastS(GtkWidget *widget, ezgl::application *app) {
    for (int numPOI = 0; numPOI < getNumPointsOfInterest(); ++numPOI) {
        float x = x_from_lon(pointofinterest[numPOI].position.lon());
        float y = y_from_lat(pointofinterest[numPOI].position.lat());
        std::string POItype = pointofinterest[numPOI].type;


        if (pointofinterest[numPOI].isFastFood == true) {
            pointofinterest[numPOI].isFastFood = false;
        } else {
            if (POItype == "fast_food") {
                pointofinterest[numPOI].isFastFood = true;
            } else {
                pointofinterest[numPOI].isFastFood = false;
            }
        }
    }
    app->refresh_drawing();
}

void walk(GtkWidget *widget, ezgl::application *app) {
    walkEnabled = true;
    GObject *window;
    GtkWidget *content_area;
    GtkWidget *label;
    GtkWidget* dialog;
    GtkWidget* entry3;
    GtkWidget* entry4;
    GtkWidget* entry5;

    window = app->get_object(app->get_main_window_id().c_str());
    dialog = gtk_dialog_new_with_buttons(("Enter Intersections"),
            (GtkWindow*) window,
            GTK_DIALOG_MODAL,
            ("Search"),
            GTK_RESPONSE_ACCEPT,
            ("Close"),
            GTK_RESPONSE_REJECT,
            NULL);
    content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
    label = gtk_label_new("Enter the names of the intersections you wish to travel between");
    gtk_container_add(GTK_CONTAINER(content_area), label);
    entry3 = gtk_entry_new();
    entry4 = gtk_entry_new();

    if (highlightCount > 0) {
        std::pair<std::string, std::string> fromTo = getNameInfo(intersections);
        gtk_entry_set_text((GtkEntry*) entry3, fromTo.first.c_str());
        if (highlightCount == 2) gtk_entry_set_text((GtkEntry*) entry4, fromTo.second.c_str());
    } else {
        gtk_entry_set_text((GtkEntry*) entry3, "From Intersection");
        gtk_entry_set_text((GtkEntry*) entry4, "To Intersection");
    }


    entry5 = gtk_entry_new();

    gtk_entry_set_text((GtkEntry*) entry5, "Walking Time (sec)-Walking speed (m/s)");
    gtk_container_add(GTK_CONTAINER(content_area), entry3);
    gtk_container_add(GTK_CONTAINER(content_area), entry4);
    gtk_container_add(GTK_CONTAINER(content_area), entry5);
    gtk_widget_show_all(dialog);

    entry1 = (GtkEntry*) entry3;
    entry2 = (GtkEntry*) entry4;
    entry6 = (GtkEntry*) entry5;


    g_signal_connect(GTK_DIALOG(dialog), "response", G_CALLBACK(search_dialog_response), app);
    app->refresh_drawing();
}
