import fiona
from fiona.crs import from_epsg
from collections import OrderedDict

output_driver = "ESRI Shapefile"
schema = {
                "geometry": "LineString",
                "properties": OrderedDict([
                    ("Flow", "str"),
                    ("Temp", "str"),
                    ("Cost", "str"),
                    ("Length", "str")
                ])
                }


def create_transmission_line_shp(transmission_lines, flows, temperatures, costs, lengths, file):
    with fiona.open(file,  "w", crs=from_epsg(4326), driver=output_driver, schema=schema) as shp:
        for transmission_line, flow, temperature, cost, length in zip(transmission_lines, flows, temperatures, costs, lengths):
            line = {
                "geometry": {
                    "type": "LineString",
                    "coordinates": transmission_line
                },
                "properties": OrderedDict([
                    ("Flow", str(flow) + " MWh/a"),
                    ("Temp", str(temperature) + " C"),
                    ("Cost", str(cost) + " Euro"),
                    ("Length", str(length) + " km")
                ])
            }
            shp.write(line)


if __name__ == "__main__":
    transmissions = [((0,0), (1,1)), ((1,1), (2,1))]
    flows = [2, 3]
    temperatures = [100, 200]
    costs = [100000, 200000]
    lengths = [4, 5]
    create_transmission_line_shp(transmissions, flows, temperatures, costs, lengths, "./test.shp")
