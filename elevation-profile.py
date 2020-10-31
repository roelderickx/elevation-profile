#!/usr/bin/env python

# elevation-profile -- generate an elevation profile for a gpx-track
#                      using a raster datasource
# Copyright (C) 2020  Roel Derickx <roel.derickx AT gmail>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import argparse, math, numpy as np
from lxml import etree
from demquery import Query
from osgeo import ogr, osr
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

# Tracks class from hikingmap, slightly modified (see MODIFIED comments)
# This class may be removed from future versions in favour of the hikingmap.Tracks class,
# but it requires porting several changes upstream
class Tracks:
    def __init__(self):
        self.tracks = list()
        self.waypoints = list()


    # Read all tracks from a given list of gpx files and store them in memory
    def parse_files(self, gpxfiles):
        for gpxfile in gpxfiles:
            print("Reading file %s" % gpxfile)

            gpxdoc = etree.parse(gpxfile)
            gpxroot = gpxdoc.getroot()

            for gpxtrack in gpxroot.findall('trk', namespaces=gpxroot.nsmap):
                self.__parse_track(gpxtrack)


    def __parse_track(self, gpxtrack):
        namenode = gpxtrack.find('name', namespaces=gpxtrack.nsmap)
        trackname = namenode.text if namenode is not None and namenode.text else "[unnamed]"
        print("Found track %s" % trackname)

        track = list()
        for coord in gpxtrack.findall('trkseg/trkpt', namespaces=gpxtrack.nsmap):
            # MODIFIED dem-query requires lon-lat tuples in stead of hikingmap.Coordinate objects
            # TODO implement to_tuple in hikingmap.Coordinate
            track.append((float(coord.get('lon')),
                          float(coord.get('lat'))))

        # search if track connects to existing track in tracks
        foundindex = 0
        foundtrack = False
        for foundindex, existingtrack in enumerate(self.tracks):
            # MODIFIED hikingmap.Coordinate.equals replaced by ==
            if existingtrack[0] == track[0]:
                print("=> same startpoint as track %d: reversing track" % foundindex)
                track.reverse()
            elif existingtrack[-1] == track[-1]:
                print("=> same endpoint as track %d: reversing track" % foundindex)
                track.reverse()

            if existingtrack[-1] == track[0]:
                print("=> connecting after track %d" % foundindex)
                newtrack = existingtrack + track[1:]
                self.tracks[foundindex] = newtrack
                foundtrack = True
                break
            elif existingtrack[0] == track[-1]:
                print("=> connecting before track %d" % foundindex)
                newtrack = track + existingtrack[1:]
                self.tracks[foundindex] = newtrack
                foundtrack = True
                break

        if not foundtrack:
            print("=> new track %d" % len(self.tracks))
            self.tracks.append(track)
    
    
    # Calculate waypoints after each waypt_distance for every track
    def calculate_waypoints(self, waypt_distance, length_unit):
        for (trackindex, track) in enumerate(self.tracks):
            # MODIFIED replaced hikingmap.Coordinate.to_string() by str()
            print("Generating waypoints for track %d: %s - %s" % \
                        (trackindex, str(track[0]), str(track[-1])))
            
            track_waypoints = list()
            cumulDistance = 0
            prev_coord = track[0]
            for coord in track:
                cumulDistance = self.__add_waypoints(track_waypoints, prev_coord, coord, \
                                                     cumulDistance, waypt_distance, length_unit)
                prev_coord = coord
            # MODIFIED add latest point
            track_waypoints.append((track[-1], "%.2f" % cumulDistance))

            print("Total track distance: %.2f %s" % (cumulDistance, length_unit))
            
            self.waypoints.append(track_waypoints)
    
    
    # MODIFIED added from hikingmap.Coordinate
    def __get_earth_radius(self, length_unit):
        if length_unit == "mi":
            return 3959
        else: # default to km
            return 6371


    # MODIFIED added from hikingmap.Coordinate
    # calculate distance in km or mi between self and coord
    def __distance_haversine(self, coord1, coord2, length_unit):
        coord1_lon = math.radians(coord1[0])
        coord1_lat = math.radians(coord1[1])
        coord2_lon = math.radians(coord2[0])
        coord2_lat = math.radians(coord2[1])
        
        dLat = coord2_lat - coord1_lat
        dLon = coord2_lon - coord1_lon

        a = math.sin(dLat/2) * math.sin(dLat/2) + \
            math.sin(dLon/2) * math.sin(dLon/2) * \
            math.cos(coord1_lat) * math.cos(coord2_lat)
        c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))

        return self.__get_earth_radius(length_unit) * c


    # MODIFIED added from hikingmap.Coordinate
    # calculate bearing between self and coord
    def __bearing(self, coord1, coord2):
        coord1_lon = math.radians(coord1[0])
        coord1_lat = math.radians(coord1[1])
        coord2_lon = math.radians(coord2[0])
        coord2_lat = math.radians(coord2[1])
        
        dLon = coord2_lon - coord1_lon

        y = math.sin(dLon) * math.cos(coord2_lat)
        x = math.cos(coord1_lat) * math.sin(coord2_lat) - \
            math.sin(coord1_lat) * math.cos(coord2_lat) * math.cos(dLon)
        return math.atan2(y, x)


    # MODIFIED added from hikingmap.Coordinate
    # returns the coordinate of the point which is on a given distance
    # from self in the direction of dest_coord
    def __calc_waypoint_on_line(self, coord1, coord2, distance, length_unit):
        coord1_lon = math.radians(coord1[0])
        coord1_lat = math.radians(coord1[1])
        coord2_lon = math.radians(coord2[0])
        coord2_lat = math.radians(coord2[1])
        
        b = self.__bearing(coord1, coord2)
        earth_radius = self.__get_earth_radius(length_unit)
        return (#lon
                math.degrees(coord1_lon + \
                             math.atan2(math.sin(b) * \
                                        math.sin(distance/earth_radius) * \
                                        math.cos(coord1_lat), \
                                        math.cos(distance/earth_radius) - \
                                        math.sin(coord1_lat) * \
                                        math.sin(coord2_lat))), \
                #lat
                math.degrees(math.asin(math.sin(coord1_lat) * \
                                       math.cos(distance/earth_radius) + \
                                       math.cos(coord1_lat) * \
                                       math.sin(distance/earth_radius) * \
                                       math.cos(b))))


    # calculate all waypoints between coord1 and coord2 and append them to track_waypoints
    # returns cumulative distance at coord2
    def __add_waypoints(self, track_waypoints, coord1, coord2, cumul_dist_at_coord1, \
                        waypt_distance, length_unit):
        # MODIFIED hikingmap.Coordinate.equals replaced by ==
        if coord1 == coord2:
            if cumul_dist_at_coord1 == 0:
                track_waypoints.append((coord1, "0.00"))
            return cumul_dist_at_coord1
        else:
            # MODIFIED hikingmap.Coordinate.distance_haversine replaced by local method
            cumul_dist_at_coord2 = \
                cumul_dist_at_coord1 + self.__distance_haversine(coord1, coord2, length_unit)
            # MODIFIED multiplication by 100 to allow for maximum 1 waypoint per 10 meter
            # TODO factor 100 is shady and slows down the algorithm
            for dist in range(int(cumul_dist_at_coord1*100) + 1, int(cumul_dist_at_coord2*100) + 1):
                if dist % (waypt_distance*100) == 0:
                    d = (dist/100) - cumul_dist_at_coord1
                    waypt = self.__calc_waypoint_on_line(coord1, coord2, d, length_unit)
                    track_waypoints.append((waypt, "%.2f" % (dist/100)))

            return cumul_dist_at_coord2



def parse_commandline():
    parser = argparse.ArgumentParser(description='Generate an elevation profile for a given GPX track')
    parser.add_argument('--datasource', dest='datasource', help='Elevation raster datafile',
                        required = True)
    parser.add_argument('--src-srs', dest='srcsrs', type=int, default=4326,
                        help='EPSG code of input data. Do not include the EPSG: prefix.')
    parser.add_argument('--gpx', dest='gpxfile',
                        help='Track for which the elevation profile should be rendered')
    params = parser.parse_args()
    
    return params


def reproject(coord, coord_trans):
    if coord_trans is None:
        return coord
    else:
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(coord[0], coord[1])
        point.Transform(coord_trans)
        return (point.GetX(), point.GetY())


def main():
    params = parse_commandline()
    
    coord_trans = None
    if params.srcsrs != 4326:
        gpx_spatial_ref = osr.SpatialReference()
        gpx_spatial_ref.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        gpx_spatial_ref.ImportFromEPSG(4326)
        
        ds_spatial_ref = osr.SpatialReference()
        ds_spatial_ref.ImportFromEPSG(params.srcsrs)
        
        coord_trans = osr.CoordinateTransformation(gpx_spatial_ref, ds_spatial_ref)

    tracks = Tracks()
    tracks.parse_files([ params.gpxfile ])
    tracks.calculate_waypoints(0.1, 'km')

    query = Query([ params.datasource ])

    for (index, track_waypoints) in enumerate(tracks.waypoints):
        print("Profile track %d" % index)
        
        # TODO generate SVG in stead
        dist_list = [ float(wpt[1]) for wpt in track_waypoints ]
        elev_list = query.query_points([ reproject(wpt[0], coord_trans) for wpt in track_waypoints ], \
                                       interp_kind='linear')
        slope_list = [ ]
        for (index, wpt) in enumerate(track_waypoints[1:]):
            slope_list.append((elev_list[index+1] - elev_list[index]) / \
                              ((dist_list[index+1] - dist_list[index]) * 10))
                              # factor 1000 for conversion from km to m, factor 100 for percentage
        slopes = np.asarray(slope_list, dtype=np.float32)
        
        # Create a set of line segments so that we can color them individually
        # This creates the points as a N x 1 x 2 array so that we can stack points
        # together easily to get the segments. The segments array for line collection
        # needs to be (numlines) x (points per line) x 2 (for x and y)
        points = np.array([dist_list, elev_list]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        
        mean_elev=round((sum(elev_list)/len(elev_list)),3)
        min_elev=min(elev_list)
        max_elev=max(elev_list)
        distance=dist_list[-1]
        
        fig, ax = plt.subplots(figsize=(50, 4), tight_layout=True)

        # Create a continuous norm to map from data points to colors
        #norm_boundary = max(abs(slopes.min()), abs(slopes.max()))
        norm = plt.Normalize(slopes.min(), slopes.max())
        colors = [ (r, 1-r, 0) for r in np.hstack((np.log(np.arange(np.e, 1, (np.e-1)/slopes.min())), \
                                                   np.log(np.arange(1, np.e, (np.e-1)/slopes.max())))) ]
        cmap = ListedColormap(colors)
        lc = LineCollection(segments, cmap=cmap, norm=norm)
        # Set the values used for colormapping
        lc.set_array(slopes)
        lc.set_linewidth(2)
        line = ax.add_collection(lc)
        fig.colorbar(line, ax=ax) # draw slope angle legend

        plt.xlabel("Distance(km)")
        plt.ylabel("Elevation(m)")
        plt.fill_between(dist_list, elev_list, 0, alpha=0.1)
        plt.xticks(np.arange(0, int(distance)+1, max(int(distance/20), 1)))
        plt.grid()

        plt.plot([0,distance], [min_elev,min_elev], '--g', label='min: %.2f m' % min_elev)
        plt.plot([0,distance], [max_elev,max_elev], '--r', label='max: %.2f m' % max_elev)
        plt.plot([0,distance], [mean_elev,mean_elev], '--y', label='ave: %.2f m' % mean_elev)
        plt.legend(fontsize='small')

        plt.show()


if __name__ == '__main__':
    main()

