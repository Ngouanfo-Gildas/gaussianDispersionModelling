from GDModeling import *
from load_data import *

print("génération de positions potentielles...")
generate_potentials_positions(100, 100, 10)
print("génération de positions potentielles terminée.")
print("génération des sources de pollution...")
generate_src_pollutions(100, 100, 1)
print("génération des sources de pollution terminée.")
print("Lecture de données...")
I = read_data_pd("pollution_source_positions.csv")
P = read_data_pd("potentials_positions.csv")
print("Lecture terminée.")
C0 = 0.000020  # (g/m^3)

print("Calcul des zones de pollution\n Patientez le temps que les caculs se font...")
Z = pollutionZone(P, I, C0)

show_genome(0, 0, Z, I)

# I_ = read_data_pd("pollution_source_positions.csv")
# P_ = read_data_pd("potentials_positions.csv")
# w_ = pollutionZone(P_, I_, C0)
# print(w_)

# lo = -0.19079100000000002
# la = 51.558191
# print(convert_long_lat_xy(long_o=0, lat_o=0, long=lo, lat=la))