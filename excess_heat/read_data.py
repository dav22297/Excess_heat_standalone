import pandas as pd
import re
import fiona
import os
import csv

from pyproj import Proj, Transformer
from shapely.geometry import Point, Polygon, MultiPolygon
from shapely.wkb import loads
import numpy as np

def extract_coordinates_from_wkb_point(point):
    """
    Function extracting the coordinates from a well known byte hexadecimal string.

    :param point: Well known byte hexadecimal string describing a point.
    :type point: string
    :return: x and y coordinate of point.
    :rtype: touple of floats.
    """
    geometry = loads(point, hex=True)
    return geometry.x, geometry.y


def ad_industrial_database_dict(dict):
    country_to_nuts0 = {"Austria": "AT", "Belgium": "BE", "Bulgaria": "BG", "Cyprus": "CY", "Czech Republic": "CZ",
                        "Germany": "DE", "Denmark": "DK", "Estonia": "EE", "Finland": "FI", "France": "FR",
                        "Greece": "EL", "Hungary": "HU", "Croatia": "HR", "Ireland": "IE", "Italy": "IT",
                        "Lithuania": "LT", "Luxembourg": "LU", "Latvia": "LV", "Malta": "MT", "Netherland": "Nl",
                        "Netherlands": "Nl",
                        "Poland": "PL", "Portugal": "PT", "Romania": "RO", "Spain": "ES", "Sweden": "SE",
                        "Slovenia": "SI", "Slovakia": "SK", "United Kingdom": "UK", "Albania": "AL", "Montenegro": "ME",
                        "North Macedonia": "MK", "Serbia": "RS", "Turkey": "TR", "Switzerland": "CH", "Iceland": "IS",
                        "Liechtenstein": "LI", "Norway": "NO"}
    raw_data = pd.DataFrame(dict["industrial_database"])
    raw_data = raw_data.loc[:, ("geom", "subsector",  "country", "excess_heat_100_200c", "excess_heat_200_500c",
                                "excess_heat_500c")]
    raw_data["Lon"] = ""
    raw_data["Lat"] = ""
    raw_data["Nuts0"] = ""
    for i, site in raw_data.iterrows():
        # check if site location is available
        if not pd.isna(site["geom"]) and site["geom"] != "":
            lon, lat = extract_coordinates_from_wkb_point(site["geom"])
            raw_data.loc[i, "Lon"] = lon
            raw_data.loc[i, "Lat"] = lat
        if not pd.isna(site["country"]) and site["country"] != "":
            raw_data.loc[i, "Nuts0"] = country_to_nuts0[site["country"]]

    raw_data = raw_data[raw_data.Lon != ""]
    raw_data = raw_data[raw_data.Lat != ""]

    data = pd.DataFrame(columns=("Lon", "Lat", "Nuts0_ID", "Subsector", "Excess_heat", "Temperature"))
    raw_data["excess_heat_100_200c"] = pd.to_numeric(raw_data["excess_heat_100_200c"])
    raw_data["excess_heat_200_500c"] = pd.to_numeric(raw_data["excess_heat_200_500c"])
    raw_data["excess_heat_500c"] = pd.to_numeric(raw_data["excess_heat_500c"])


    for i, site in raw_data.iterrows():
        # check if heat at specific temperature range is available
        # TODO deal with units; hard coded temp ranges?
        if not pd.isna(site["excess_heat_100_200c"]) and site["excess_heat_100_200c"] != "" and site["excess_heat_100_200c"] != 0:
            data.loc[data.shape[0]] = (site["Lon"], site["Lat"], site["Nuts0"], site["subsector"],
                                       1000*site["excess_heat_100_200c"], 150)
        if not pd.isna(site["excess_heat_200_500c"]) and site["excess_heat_200_500c"] != "" and site["excess_heat_200_500c"] != 0:
            data.loc[data.shape[0]] = (site["Lon"], site["Lat"], site["Nuts0"],
                                       site["subsector"], 1000*site["excess_heat_200_500c"], 350)
        if not pd.isna(site["excess_heat_500c"]) and site["excess_heat_500c"] != "" and site["excess_heat_500c"] != 0:
            data.loc[data.shape[0]] = (site["Lon"], site["Lat"], site["Nuts0"],
                                       site["subsector"], 1000*site["excess_heat_500c"], 500)

    return data


def ad_TUW23(out_shp_label, nuts2_id):
    """
    Function extracting potential heat sinks computed by the TUW23 CM. It creates a grid of points of constant density
    inside coherent areas.

    :param out_shp_label: File name of shp file containing the coherent areas of TUW23 CM.
    :type out_shp_label: sting
    :return: Dataframe containing the potential heat sinks and a correspondence id for each coherent aera.
    :rtype: pandas Dataframe
    """
    try:
        coherent_areas = fiona.open(out_shp_label)
    except IOError:
        return -1
    inProj = Proj(init='epsg:3035')
    outProj = Proj(init='epsg:4326')
    transformer = Transformer.from_proj(inProj, outProj)

    coherent_areas_transformed = []
    for coherent_area in coherent_areas:
        coordinates = []
        if coherent_area["geometry"]["type"] == "Polygon":
            for coordinate in coherent_area["geometry"]["coordinates"][0]:
                coordinates.append(transformer.transform(*coordinate))
            poly = Polygon(coordinates)
            poly.heat_dem = 1000*float(re.findall("\d+\.\d+", coherent_area["properties"]["Potential"])[0])
            coherent_areas_transformed.append(poly)
        elif coherent_area["geometry"]["type"] == "MultiPolygon":
            multipolygon_polygons = []
            for polygon in coherent_area["geometry"]["coordinates"][0]:
                coordinates = []
                for coordinate in polygon:
                    coordinates.append(transformer.transform(*coordinate))
                multipolygon_polygons.append(Polygon(coordinates))
            multi_poly = MultiPolygon(multipolygon_polygons)
            multi_poly.heat_dem = 1000*float(re.findall("\d+\.\d+", coherent_area["properties"]["Potential"])[0])
            coherent_areas_transformed.append(multi_poly)

    data = []
    delta = 0.015
    for i, coherent_area in enumerate(coherent_areas_transformed):
        entry_points = []
        (minx, miny, maxx, maxy) = coherent_area.bounds
        possible_x = np.arange(minx, maxx, delta)
        possible_y = np.arange(miny, maxy, delta)
        for x in possible_x:
            for y in possible_y:
                if coherent_area.buffer(0).contains(Point(x, y)):
                    entry_points.append((x, y))
        if len(entry_points) == 0:
            entry_points.append((maxx, maxy))
        induvidual_heat_demand = coherent_area.heat_dem / len(entry_points)
        for entry_point in entry_points:
            data.append([*entry_point, induvidual_heat_demand, -i])

    data = pd.DataFrame(data, columns=["Lon", "Lat", "Heat_demand", "id"])
    data["Nuts2_ID"] = nuts2_id

    data["ellipsoid"] = "SRID=4326"
    data["Economic_Activity"] = "Steam and air conditioning supply"
    data["Temperature"] = 100

    return data


def ad_industry_profiles_dict(dicts):
    dict_names = ["load_profile_industry_chemicals_and_petrochemicals_yearlong_2018", "load_profile_industry_food_and_tobacco_yearlong_2018",
                  "load_profile_industry_iron_and_steel_yearlong_2018", "load_profile_industry_non_metalic_minerals_yearlong_2018",
                  "load_profile_industry_paper_yearlong_2018"]
    data = []
    for name, dict in zip(dict_names, dicts):
        raw_data = pd.DataFrame(dict[name])
        raw_data = raw_data.loc[:, ("NUTS0_code", "process", "hour", "load")]
        raw_data["load"] = pd.to_numeric(raw_data["load"])
        raw_data["hour"] = pd.to_numeric(raw_data["hour"])
        data.append(raw_data)

    return data


def ad_residential_heating_profile_dict(dict):

    data = pd.DataFrame(dict["load_profile_residential_heating_yearlong_2010"])
    data = data.loc[:, ("NUTS2_code", "process", "hour", "load")]
    data["load"] = pd.to_numeric(data["load"])
    data["hour"] = pd.to_numeric(data["hour"])
    return data


def ad_industry_profiles_local(nuts0_ids):
    """
    Loads industry profiles of different subcategories from different csv files.

    :return: List of dataframes containing the csv files data.
    :rtype: list [pd.Dataframe, pd.Dataframe, ...].
    """

    file_names = ("hotmaps_task_2.7_load_profile_industry_chemicals_and_petrochemicals_yearlong_2018.csv",
                  "hotmaps_task_2.7_load_profile_industry_food_and_tobacco_yearlong_2018.csv",
                  "hotmaps_task_2.7_load_profile_industry_iron_and_steel_yearlong_2018.csv",
                  "hotmaps_task_2.7_load_profile_industry_non_metalic_minerals_yearlong_2018.csv",
                  "hotmaps_task_2.7_load_profile_industry_paper_yearlong_2018.csv")

    path = os.path.dirname(
           os.path.dirname(os.path.abspath(__file__)))
    path = os.path.join(path, "data")

    data = []
    for file_name in file_names:
        sub_path = os.path.join(path, file_name)
        # determine delimiter of csv file
        with open(sub_path, 'r', encoding='utf-8') as csv_file:
            delimiter = csv.Sniffer().sniff(csv_file.readline()).delimiter

        raw_data = pd.read_csv(sub_path, sep=delimiter, usecols=("NUTS0_code", "process", "hour", "load"))
        raw_data = raw_data[raw_data["NUTS0_code"].isin(nuts0_ids)]
        data.append(raw_data)

    return data


def ad_residential_heating_profile_local(nuts2_ids):
    """
    Loads residential heating profiles from csv file.

    :return: Dataframe containing the data of the csv file.
    :rtype: pandas dataframe.
    """

    path = os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))
    path = os.path.join(path, "data")
    path1 = os.path.join(path, "hotmaps_task_2.7_load_profile_residential_heating_yearlong_2010_part1.csv")
    path2 = os.path.join(path, "hotmaps_task_2.7_load_profile_residential_heating_yearlong_2010_part2.csv")

    # determine delimiter of csv file
    with open(path1, 'r', encoding='utf-8') as csv_file:
        delimiter = csv.Sniffer().sniff(csv_file.readline()).delimiter
    data = pd.read_csv(path1, sep=delimiter, usecols=("NUTS2_code", "process", "hour", "load"))
    with open(path2, 'r', encoding='utf-8') as csv_file:
        delimiter = csv.Sniffer().sniff(csv_file.readline()).delimiter
    data2 = pd.read_csv(path2, sep=delimiter, usecols=("NUTS2_code", "process", "hour", "load"))

    data = data.append(data2)
    data = data[data["NUTS2_code"].isin(nuts2_ids)]

    return data


def ad_industrial_database_local(nuts0_ids):
    """
    loads data of heat sources given by a csv file.

    :return: dataframe containing the data of the csv file.
    :rtype: pandas dataframe.
    """

    country_to_nuts0 = {"Austria": "AT", "Belgium": "BE", "Bulgaria": "BG", "Cyprus": "CY", "Czech Republic": "CZ",
                    "Germany": "DE", "Denmark": "DK", "Estonia": "EE", "Finland": "FI", "France": "FR",
                    "Greece": "EL", "Hungary": "HU", "Croatia": "HR", "Ireland": "IE", "Italy": "IT",
                    "Lithuania": "LT", "Luxembourg": "LU", "Latvia": "LV", "Malta": "MT", "Netherland": "Nl",
                    "Netherlands": "Nl",
                    "Poland": "PL", "Portugal": "PT", "Romania": "RO", "Spain": "ES", "Sweden": "SE",
                    "Slovenia": "SI", "Slovakia": "SK", "United Kingdom": "UK", "Albania": "AL", "Montenegro": "ME",
                    "North Macedonia": "MK", "Serbia": "RS", "Turkey": "TR", "Switzerland": "CH", "Iceland": "IS",
                    "Liechtenstein": "LI", "Norway": "NO"}
    path = os.path.dirname(
           os.path.dirname(os.path.abspath(__file__)))
    path = os.path.join(os.path.join(path, "data"), "Industrial_Database.csv")

    # determine delimiter of csv file
    with open(path, 'r', encoding='utf-8') as csv_file:
        delimiter = csv.Sniffer().sniff(csv_file.readline()).delimiter

    raw_data = pd.read_csv(path, sep=delimiter, usecols=("geom", "Subsector", "Excess_Heat_100-200C",
                                                         "Excess_Heat_200-500C", "Excess_Heat_500C", "Country"))

    # dataframe for processed data
    data = pd.DataFrame(columns=("ellipsoid", "Lon", "Lat", "Nuts0_ID", "Subsector", "Excess_heat", "Temperature"))
    for i, site in raw_data.iterrows():
        # check if site location is available
        if not pd.isna(site["geom"]):
            # extract ellipsoid model and (lon, lat) from the "geom" column
            ellipsoid, coordinate = site["geom"].split(";")
            m = re.search("[-+]?[0-9]*\.?[0-9]+.[-+]?[0-9]*\.?[0-9]+", coordinate)
            m = m.group(0)
            lon, lat = m.split(" ")
            lon = float(lon)
            lat = float(lat)

            nuts0 = country_to_nuts0[site["Country"]]

            # check if heat at specific temperature range is available
            # TODO deal with units; hard coded temp ranges?
            if not pd.isna(site["Excess_Heat_100-200C"]) and site["Excess_Heat_100-200C"] != "" and site["Excess_Heat_100-200C"] != 0:
                data.loc[data.shape[0]] = (ellipsoid, lon, lat, nuts0, site["Subsector"],
                                           site["Excess_Heat_100-200C"] * 1000, 150)
            if not pd.isna(site["Excess_Heat_200-500C"]) and site["Excess_Heat_200-500C"] != "" and site["Excess_Heat_200-500C"] != 0:
                data.loc[data.shape[0]] = (ellipsoid, lon, lat, nuts0,
                                           site["Subsector"], site["Excess_Heat_200-500C"] * 1000, 350)
            if not pd.isna(site["Excess_Heat_500C"]) and site["Excess_Heat_500C"] != "" and site["Excess_Heat_500C"] != 0:
                data.loc[data.shape[0]] = (ellipsoid, lon, lat, nuts0,
                                           site["Subsector"], site["Excess_Heat_500C"] * 1000, 500)

    data = data[data["Nuts0_ID"].isin(nuts0_ids)]

    return data


def ad_entry_points(nuts2_id):
    """
    loads data of heat sources given by a csv file.

    :return: dataframe containing the data of the csv file.
    :rtype: pandas dataframe.
    """

    path = os.path.dirname(
           os.path.dirname(os.path.abspath(__file__)))
    path = os.path.join(os.path.join(path, "data"), "entry_points.csv")

    # determine delimiter of csv file
    with open(path, 'r', encoding='utf-8') as csv_file:
        delimiter = csv.Sniffer().sniff(csv_file.readline()).delimiter
    raw_data = pd.read_csv(path, sep=delimiter, usecols=("Lon", "Lat", "Annual heat demand in Gwh", "id"))
    data = pd.DataFrame(raw_data, columns=("Lon", "Lat", "Heat_demand", "id"))
    data["Heat_demand"] = 1000 * data["Heat_demand"]
    data["Nuts2_ID"] = nuts2_id

    return data
