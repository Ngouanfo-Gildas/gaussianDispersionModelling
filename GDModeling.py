import math
import numpy as np

""" Gaussian Dispersion Modeling implemention"""

# paramètres
Q = 5         # (g/s) taux d'émission de la source de pollution
hs = 25       # (m) hauteur de sources de pollution
Vw = 5        # (m/s) vitesse du vent
Dw = 0      # direction du vent (°)
Ts = 303.15     # (30°c = ?K) température du polluant à la source
T = 280.15      # (7°c = ?K) température de l'air ambiante
g = 9.8         # (m/s²) constante de gravité
ay = 1.36
by = 0.82
az = 0.275
bz = 0.69
V = 0.0000000019  # (m³/s)

# Calcul de F; 
# ON SUPPOSE QUE Ts > T
F = (g/math.pi)*V*((Ts-T)/Ts)  # = 0.31350340690241546

def delta_h(F, Vw, x):  # (en mètre)
    f = pow(F, 1/3.0)
    s = pow(abs(x), 2/3.0)
    l = 1.6 * f * s
    return l/Vw

def calcsigmas(x1):

    x=np.abs(x1)

    a=np.zeros(np.shape([x]))
    b=np.zeros(np.shape([x]))
    c=np.zeros(np.shape([x]))
    d=np.zeros(np.shape([x]))

    # vertical
    ind=np.where((x<100.) & (x>0.))
    a=122.800
    b=0.94470
    
    ind=np.where((x>=100.) & (x<150.))
    a=158.080
    b=1.05420
		
    ind=np.where((x>=150.) & (x<200.))
    a=170.220
    b=1.09320
    
    ind=np.where((x>=200.) & (x<250.))
    a=179.520
    b=1.12620
		
    ind=np.where((x>=250.) & (x<300.))
    a=217.410
    b=1.26440
    
    ind=np.where((x>=300.) & (x<400.))
    a=258.89
    b=1.40940
		
    ind=np.where((x>=400.) & (x<500.))
    a=346.75
    b=1.7283
    
    ind=np.where((x>=500.) & (x<3110.))
    a=453.85
    b=2.1166
    
    ind=np.where((x>=3110.))
    a=453.85
    b=2.1166
		
    # cross wind
    c[:]=24.1670
    d[:]=2.5334

    sig_z = a*(x/1000.)**b
    sig_z[np.where(sig_z[:]>5000.)] = 5000.

    theta = 0.017453293*(c-d*np.log(np.abs(x+1e-15)/1000.))
    sig_y = 465.11628*x/1000.*np.tan(theta)
    
    return sig_y, sig_z

# ON SUPPOSE QUE x > 0
# coeficient de dispersion horizontal et vertical
sigma_y = lambda x, ay, by : ay * pow(abs(x), by)
sigma_z = lambda x, az, bz : az * pow(abs(x), bz)

"""print("\nsig_y pm = ", sigma_y(10, ay, by), "\nsig_z pm = ",sigma_z(10, az, bz))
sig_y, sig_z = calcsigmas(10)
print("\nsig_y = ", sig_y, "\nsig_z = ",sig_z)"""


def gaussianDM(x, y, z, Vw, Q, hs, ay, az, by, bz):
    """ Concentration du polluant au point (0, 0, hs) """
    Delta_h = delta_h(F, Vw, x)
    H = hs + Delta_h
    z = hs
    sigmaY = sigma_y(x, ay, by)
    sigmaZ = sigma_z(x, az, bz) 

    # sigmaY, sigmaZ = calcsigmas(10)

    """ C = Q * v1* v2 / v3 """
    sigY2 = 2 * sigmaY       #sigmaY²
    v1 = math.exp(-pow(y, 2) / sigY2)
    v2 = math.exp(-pow((z-H), 2) / sigY2) +\
        math.exp(-pow((z+H), 2) / sigY2)
    v3 = math.pi * Vw * sigY2 * sigmaZ
    # print("x, y, sigmaY, sigmaZ, v1, v2, v3",x , y, sigmaY, sigmaZ, v1, v2, v3)
    d1 = v1 * v2
    d2 = Q / v3
    return  d1 * d2


def gauss_model_func(x, y, z, Vw, Q, xs, ys, hs, wind_dir):
    """wind_dir = direction du vent en dégré\n
        u = vitesse du vent en m/s\n
        x,y,z = position de mesure\n
        xs,ys,H = source de polution
        """
    # Dy = 10.
    # Dz = 10.
    # STABILITY = 4
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

    # ind = np.where(downwind>0.)
    # C = np.zeros(len(x),len(y))
    
    # calculate sigmas based on stability and distance downwind
    # (sig_y,sig_z) = calcsigmas(downwind)
    sig_y = sigma_y(downwind, ay, by)
    sig_z = sigma_z(downwind, az, bz)

    Delta_h = delta_h(F, Vw, downwind)
    H = hs + Delta_h

    # C = Q/(2.*np.pi*u1*sig_y*sig_z) \
    #     * np.exp(-crosswind**2./(2.*sig_y**2.))  \
    #     *(np.exp(-(z-H)**2./(2.*sig_z**2.)) + \
    #     np.exp(-(z+H)**2./(2.*sig_z**2.)))
    C = Q/(2.*np.pi*Vw*sig_y*sig_z) \
        * np.exp(-crosswind**2./(2.*sig_y)) \
        * (np.exp(-(z-H)**2./(2.*sig_y)) \
        + np.exp(-(z+H)**2./(2.*sig_y)))
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
    M = len(I)
    N = len(P)
    W = np.zeros((M,N), dtype=int)
    Z = list()
    an = 0
    for i in range(M):
        an += 1
        xs = I[i][0]
        ys = I[i][1]
        for p in range(N):
            x = P[p][0]
            y = P[p][1]
            # x, y = convert_long_lat_xy(long_o=0, lat_o=0, long=P[p][0], lat=P[p][1])
            z = 20  # hauteur du point potentiel
            C = gauss_model_func(x, y, z, Vw, Q, xs, ys, hs, Dw)
            if C >= C0:
                W[i][p] = 1
                Z.append(P[p])
                print("Un nouveau point trouvé: ", P[p], "pour la source", I[i])
                # Z.append((I[i], P[p]))
    return Z
