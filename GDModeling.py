import math
import numpy as np

""" Gaussian Dispersion Modeling implemention"""

### Paramètres du modèle Gaussien
hs = 25           # (m) hauteur de sources de pollution
Q = 5             # (g/s) taux d'émission de la source de pollution
V = 0.0000000019  # (m³/s)
Ts = 303.15       # (30°c = ?K) température du polluant à la source
ay = 1.36         
by = 0.82 
az = 0.275 
bz = 0.69 
g = 9.8           # (m/s²) constante de gravité
### Paramètres météorologiques
Vw = 5      # (m/s) vitesse du vent
Dw = 35     # direction du vent (°)
T = 280.15  # (7°c = ?K) température de l'air ambiante

# Calcul de F; ON SUPPOSE QUE Ts > T
F = (g/math.pi)*V*((Ts-T)/Ts)  
# = 0.31350340690241546

def delta_h(F, Vw, x):  # (en mètre)
    f = pow(F, 1/3.0)
    s = pow(abs(x), 2/3.0)
    l = 1.6 * f * s
    return l/Vw
    
# ON SUPPOSE QUE x > 0
# coeficient de dispersion horizontal et vertical
sigma_y = lambda x, ay, by : ay * pow(abs(x), by)
sigma_z = lambda x, az, bz : az * pow(abs(x), bz)

def gauss_model_func(x, y, z, Vw, Q, xs, ys, hs, wind_dir):
    """wind_dir = direction du vent en dégré\n
        Vw = vitesse du vent en m/s\n
        x,y,z = position de mesure\n
        xs,ys,H = source de polution
        """
    u1 = Vw
    # déplacer les coordonnées pour que la pile soit le point central
    x1 = x - xs; # shift the coordinates so that stack is centre point
    y1 = y - ys; 

    # components of u in x and y directions (direction du vent)
    wx = u1 * np.sin((wind_dir - 180.) * np.pi / 180.)
    wy = u1 * np.cos((wind_dir - 180.) * np.pi / 180.)

    # Need angle between point x, y and the wind direction, so use scalar product:
    dot_product = wx * x1 + wy * y1
    # product of magnitude of vectors:
    magnitudes = u1 * np.sqrt(x1**2. + y1**2.)

    # angle between wind and point (x,y)
    subtended = np.arccos(dot_product / (magnitudes + 1e-15))
    # distance to point x,y from stack
    hypotenuse = np.sqrt(x1**2. + y1**2.)

    # distance along the wind direction to perpendilcular line that intesects
    # x,y
    downwind = np.cos(subtended) * hypotenuse

    # Now calculate distance cross wind.
    crosswind = np.sin(subtended) * hypotenuse

    sig_y = sigma_y(downwind, ay, by)
    sig_z = sigma_z(downwind, az, bz)

    Delta_h = delta_h(F, Vw, downwind)
    H = hs + Delta_h

    C = Q/(np.pi*Vw*2*sig_y*sig_z) * np.exp(-crosswind**2/(2*sig_y)) \
        * (np.exp(-(z-H)**2/(2*sig_y)) + np.exp(-(z+H)**2/(2*sig_y)))
    return C

def convert_long_lat2xy(long_o, lat_o, long, lat):
    """
        Converti (longitude, latitude) en coordonnées (X, Y) du plan
    """
    dx = (long_o - long) * 40000 * math.cos((lat_o - lat) * math.pi / 360) / 360
    dy = (lat_o - lat) * 40000 / 360
    return dx, dy

# Calcul des Zi
def pollutionZone(P, I, C0):
    M = len(I)   # sources de pollution
    N = len(P)   # positions potentielles
    W = np.zeros((M,N), dtype=int)
    Zi = list()
    Zones = list()
    z = 20  # hauteur du point potentiel
    an = 0
    result_file = open("./result_file.txt", "a")
    for i in range(M):
        H = list()
        an += 1
        xs = I[i][0]
        ys = I[i][1]
        for p in range(N):
            x = P[p][0]
            y = P[p][1]
            # x, y = convert_long_lat_xy(long_o=0, lat_o=0, long=P[p][0], lat=P[p][1])
            C = gauss_model_func(x, y, z, Vw, Q, xs, ys, hs, Dw)
            if C >= C0:
                W[i][p] = 1
                Zi.append(P[p])
                H.append((P[p], C))
        Zones.append(H)
        result_file.write(str(H)+"\n")
        print("point {}: {}".format(i ,H))
    result_file.close()
    return Zi 
