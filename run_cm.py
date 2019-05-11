from excess_heat.excess_heat import excess_heat

###########################################
# Modify by user
search_radius = 20  # km
investment_period = 20  # years
transmission_line_threshold = 0.5  # ct/kWh/a
nuts2_id = "DK05"
###########################################


district_heating_shp_file = "./data/district_heating_shp.shp"
output = "./results/results"

excess_heat(district_heating_shp_file, search_radius, investment_period, transmission_line_threshold, nuts2_id, output)
