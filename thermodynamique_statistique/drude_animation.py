#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#GlowScript 3.0 VPython

# Hard-sphere gas.

# Bruce Sherwood
# Claudine Allen
"""

from vpython import *
import numpy as np
from scipy.stats import maxwell
import itertools
from copy import deepcopy
from scipy.constants import m_e

# win = 500 # peut aider à définir la taille d'un autre objet visuel comme un histogramme proportionnellement à la taille du canevas.

# Déclaration de variables influençant le temps d'exécution de la simulation
Natoms = 50  # change this to have more or fewer atoms
dt = 1E-7  # pas d'incrémentation temporel

# Déclaration de variables physiques "Typical values"
DIM = 2 #Nombre de degrés de liberté de la simulation 
mass = m_e
Ratom = 0.01 # wildly exaggerated size of an atom
k = 1.4E-23 # Boltzmann constant
T = 10 # around room temperature
E = vector(0,0,0) * 1e-22 # electric field

#### CANEVAS DE FOND ####
L = 1 # container is a cube L on a side
gray = color.gray(0.7) # color of edges of container and spheres below
animation = canvas( width=750, height=500) # , align='left')
animation.range = L
# animation.title = 'Théorie cinétique des gaz parfaits'
# s = """  Simulation de particules modélisées en sphères dures pour représenter leur trajectoire ballistique avec collisions. Une sphère est colorée et grossie seulement pour l’effet visuel permettant de suivre sa trajectoire plus facilement dans l'animation, sa cinétique est identique à toutes les autres particules.

# """
# animation.caption = s

#### ARÊTES DE BOÎTE 2D ####
d = L/2+Ratom
r = 0.005
cadre = curve(color=gray, radius=r)
cadre.append([vector(-d,-d,0), vector(d,-d,0), vector(d,d,0), vector(-d,d,0), vector(-d,-d,0)])

#### POSITION ET QUANTITÉ DE MOUVEMENT INITIALE DES SPHÈRES ####
Atoms = [] # Objet qui contiendra les sphères pour l'animation
p = [] # quantité de mouvement des sphères
apos = [] # position des sphères
pavg = sqrt(2*mass*(DIM/2)*k*T) # average kinetic energy in 3D p**2/(2mass) = (3/2)kT : Principe de l'équipartition de l'énergie en thermodynamique statistique classique
complete_p_mag = []
complete_apos = []

ionic_cores_rows, ionic_cores_cols = 5, 5
ionic_cores_x, ionic_cores_y = np.meshgrid(
    np.linspace(-d, d, ionic_cores_cols+2)[1:-1],
    np.linspace(-d, d, ionic_cores_rows+2)[1:-1]
)   # +2 to remove cores that would be on the edges of the box
ionic_cores_func = np.vectorize(lambda x, y: simple_sphere(pos=vector(x,y,0), radius=Ratom*2, color=vec(0.82,0.72,1)))
ionic_cores = ionic_cores_func(ionic_cores_x, ionic_cores_y).flatten()

for i in range(Natoms):
    x = L*random()-L/2 # position aléatoire qui tient compte que l'origine est au centre de la boîte
    y = L*random()-L/2
    z = 0
    if i == 0:  # garde une sphère plus grosse et colorée parmis toutes les grises
        Atoms.append(simple_sphere(pos=vector(x,y,z), radius=0.03, color=color.magenta)) #, make_trail=True, retain=100, trail_radius=0.3*Ratom))
    else: Atoms.append(simple_sphere(pos=vector(x,y,z), radius=Ratom, color=gray))
    apos.append(vec(x,y,z)) # liste de la position initiale de toutes les sphères
#    theta = pi*random() # direction de coordonnées sphériques, superflue en 2D
    phi = 2*pi*random() # direction aléatoire pour la quantité de mouvement
    px = pavg*cos(phi)  # qte de mvt initiale selon l'équipartition
    py = pavg*sin(phi)
    pz = 0
    p.append(vector(px,py,pz)) # liste de la quantité de mouvement initiale de toutes les sphères

def check_ionic_core_collisions():
    hitlist = []
    r2 = 2*Ratom
    for atom_i, ionic_core in itertools.product(range(len(Atoms)), ionic_cores):
        if (apos[atom_i] - ionic_core.pos).mag < r2:
            hitlist.append(atom_i)
    return hitlist

mb_distrib_scale = np.sqrt(k * T / mass)

def apply_electric_field(hitlist: list):
    global p
    for i in range(len(p)):
        if i not in hitlist:
            p[i] += - E * dt

#### BOUCLE PRINCIPALE POUR L'ÉVOLUTION TEMPORELLE DE PAS dt ####
## ATTENTION : la boucle laisse aller l'animation aussi longtemps que souhaité, assurez-vous de savoir comment interrompre vous-même correctement (souvent `ctrl+c`, mais peut varier)
## ALTERNATIVE : vous pouvez bien sûr remplacer la boucle "while" par une boucle "for" avec un nombre d'itérations suffisant pour obtenir une bonne distribution statistique à l'équilibre
for i in range(10000):
    # rate(300)  # limite la vitesse de calcul de la simulation pour que l'animation soit visible à l'oeil humain!

    #### DÉPLACE TOUTES LES SPHÈRES D'UN PAS SPATIAL deltax
    vitesse = []   # vitesse instantanée de chaque sphère
    deltax = []  # pas de position de chaque sphère correspondant à l'incrément de temps dt
    for i in range(Natoms):
        vitesse.append(p[i]/mass)   # par définition de la quantité de nouvement pour chaque sphère
        deltax.append(vitesse[i] * dt)   # différence avant pour calculer l'incrément de position
        Atoms[i].pos = apos[i] = apos[i] + deltax[i]  # nouvelle position de l'atome après l'incrément de temps dt

    #### CONSERVE LA QUANTITÉ DE MOUVEMENT AUX COLLISIONS AVEC LES MURS DE LA BOÎTE ####
    for i in range(Natoms):
        loc = apos[i]
        if abs(loc.x) > L/2:
            if loc.x < 0: p[i].x =  abs(p[i].x)  # renverse composante x au mur de gauche
            else: p[i].x =  -abs(p[i].x)   # renverse composante x au mur de droite
        if abs(loc.y) > L/2:
            if loc.y < 0: p[i].y = abs(p[i].y)  # renverse composante y au mur du bas
            else: p[i].y =  -abs(p[i].y)  # renverse composante y au mur du haut

    ### LET'S FIND THESE COLLISIONS!!! ####
    hitlist = check_ionic_core_collisions()

    #### CONSERVE LA QUANTITÉ DE MOUVEMENT AUX COLLISIONS ENTRE SPHÈRES ####
    for i in hitlist:
        p[i] = vector(maxwell.rvs(mb_distrib_scale) * mass, 0, 0).rotate(angle=np.random.random() * 2*np.pi)

    apply_electric_field(hitlist)

    complete_p_mag.append(deepcopy(p))
    complete_apos.append(deepcopy(apos))


# Save necessary data
# with open("thermodynamique_statistique/data/drude_p.csv", "w") as f:
#     for p_step in complete_p_mag:
#         f.write(",".join(map(lambda p: f"{p.x:.4e},{p.y:.4e}", p_step)) + "\n")
# with open(f"thermodynamique_statistique/data/drude_apos_E={E.y:.0e}.csv", "w") as f:
#     for x_step in complete_apos:
#         f.write(",".join(map(lambda v: f"{v.x:.4f},{v.y:.4f}", x_step)) + "\n")
