import numpy as np
from numpy import linalg as LA
from random import random
import itertools


Natoms = 200
dt = 1E-5
n_dim = 2
mass = 4E-3/6E23 # helium mass
Ratom = 0.01
k_B = 1.4E-23 # Boltzmann constant
T = 300 # around room temperature
L = 1 # container is a cube L on a side

avg_p = np.sqrt(2*mass*(n_dim/2)*k_B*T) # average kinetic energy in 3D p**2/(2mass) = (3/2)kT : Principe de l'équipartition de l'énergie en thermodynamique statistique classique
phi = np.random.rand(Natoms)
p = avg_p * np.stack([np.cos(phi), np.sin(phi), 0*phi], 1) # quantité de mouvement des sphères
pos = np.stack([np.array([L*random()-L/2, L*random()-L/2, 0]) for _ in range(Natoms)], 0) # position des sphères

target_particle_index = 0 # index of the particle whose displacement will be saved at each colision
targeted_particle = np.array([[0.0, 0.0, 0.0, 0.0]]) # starting array to save the targeted particle's displacement
def check_collisions():
    hitlist = []
    Datom = 2*Ratom
    for i, j in itertools.combinations([i for i in range(Natoms)], 2): # Checks non-duplicate pairs of particles
        if (pos[i, 0] - pos[j,0])**2+(pos[i,1] - pos[j,1])**2+(pos[i,2] - pos[j,2])**2 < Datom**2: # 
        # if LA.norm(pos[i]-pos[j]) < Datom:
            hitlist.append((i,j))
    return hitlist

def follow_particle(data, particle_index, step, time_step, end_path=False):
    if end_path:
        return data[2:, :]
    collision_data = np.concatenate(
        (np.array([step*time_step-data[-1,0]]), p[particle_index]/mass*(step*time_step-data[-1,0])),
        axis=0
    )
    return np.vstack((data, collision_data))

for step in range(10000):
    pos += p*dt/mass

    for i in range(Natoms):
        if  abs(pos[i][0]) > L/2:
            p[i][0] *= -1
        if  abs(pos[i][1]) > L/2:
            p[i][1] *= -1
    hitlist = check_collisions()
    prev_particles = []
    #### CONSERVE LA QUANTITÉ DE MOUVEMENT AUX COLLISIONS ENTRE SPHÈRES ####
    for i, j in hitlist:
        if target_particle_index in (i,j):
            if target_particle_index not in prev_particles:
                targeted_particle = follow_particle(targeted_particle, target_particle_index, step, dt) # To compute the inter-collision speed, a maximum of one collision per time step should be used
            prev_particles += [i, j]
        # defining new variables for each pair of particles interacting
        vcm = (p[i]+p[j])/(2*mass)
        posi, posj = pos[i], pos[j]
        vi, vj = p[i]/mass, p[j]/mass
        rrel, vrel = posi-posj, vj-vi  # vecteur pour la distance entre les centres des 2 sphères
        # exclusion de cas où il n'y a pas de changements à faire
        if LA.norm(vrel) == 0: continue  # exactly same velocities si et seulement si le vecteur vrel devient nul, la trajectoire des 2 sphères continue alors côte à côte
        if LA.norm(rrel) > Ratom: continue  # one atom went all the way through another, la collision a été "manquée" à l'intérieur du pas deltax

        # calcule la distance et temps d'interpénétration des sphères dures qui ne doit pas se produire dans ce modèle
        dx = np.dot(rrel, vrel)/LA.norm(vrel)        # length of the component of rrel parallel to vrel
        dy = LA.norm(np.cross(rrel, vrel)/LA.norm(vrel)) # length of the component of rrel perpendicular to vrel
        alpha = np.asin(dy/(2*Ratom))       # alpha is the angle of the triangle composed of rrel, path of atom j, and a line from the center of atom i to the center of atom j where atome j hits atom i
        d = (2*Ratom)*np.cos(alpha)-dx      # distance traveled into the atom from first contact
        deltat = d/LA.norm(vrel) # time spent moving from first contact to position inside atom
        #### CHANGE L'INTERPÉNÉTRATION DES SPHÈRES PAR LA CINÉTIQUE DE COLLISION ####
        posi_1, posj_1 = posi-vi*deltat, posj-vj*deltat # back up to contact configuration
        pcmi, pcmj = p[i]-mass*vcm, p[j]-mass*vcm   # transform momenta to center-of-momentum frame
        pcmi, pcmj = pcmi-2*rrel*np.dot(pcmi, rrel)/LA.norm(rrel)**2, pcmj-2*rrel*np.dot(pcmj, rrel)/(LA.norm(rrel)**2) # bounce in center-of-momentum frame
        p[i], p[j] = pcmi+mass*vcm, pcmj+mass*vcm # transform momenta back to lab frame
        pos[i], pos[j] = posi+(p[i]/mass)*deltat, posj+(p[j]/mass)*deltat # move forward deltat in time, ramenant au même temps où sont rendues les autres sphères dans l'itération
targeted_particle = follow_particle(targeted_particle, target_particle_index, step, dt, end_path=True)
