import numpy as np
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
p = np.stack([avg_p*np.array([np.cos(2*np.pi*random()), np.sin(2*np.pi*random()), 0]) for _ in range(Natoms)], 0) # quantité de mouvement des sphères
pos = np.stack([np.array([L*random()-L/2, L*random()-L/2, 0]) for _ in range(Natoms)], 0) # position des sphères

target_particle_index = 0
targeted_particle = np.array([[0.0, 0.0, 0.0, 0.0]])

def check_collisions():
    hitlist = []
    for i, j in itertools.combinations([i for i in range(Natoms)], 2):
        if np.sum((pos[i] - pos[j])**2) < (2*Ratom)**2:
            hitlist.append((i,j))
    return hitlist

def follow_particle(data, particle_index, step, time_step, clean=False):
    if clean:
        return data[:-1, :2]
    time = step*time_step
    if np.all(np.isclose(data[-1,:], np.array([[0.0, 0.0, 0.0, 0.0]]), time_step) == True):
        data =  np.expand_dims(np.concatenate((np.array([time]), pos[particle_index])), axis=0)
    else:
        data[-1, :] = np.array([time-data[-1, 0], np.sqrt(np.sum((pos[particle_index] - data[-1, 1:])**2)), 0, 0])
        data =  np.concatenate(
            (data, np.expand_dims(np.concatenate((np.array([time]), pos[particle_index])), axis=0)),
            0
        )
    return data

for step in range(100):
    pos = [pos[i] + p[i]*dt/mass for i in range(Natoms)]

    for i in range(Natoms):
        if  abs(pos[i][0]) > L/2:
            p[i][0] *= -1
        if  abs(pos[i][1]) > L/2:
            p[i][1] *= -1

    hitlist = check_collisions()

    #### CONSERVE LA QUANTITÉ DE MOUVEMENT AUX COLLISIONS ENTRE SPHÈRES ####
    for i, j in hitlist:
        if target_particle_index in (i,j):
            targeted_particle = follow_particle(targeted_particle, target_particle_index, step, dt)
        # définition de nouvelles variables pour chaque paire de sphères en collision
        vcm = (p[i]+p[j])/(2*mass)
        posi, posj = pos[i], pos[j]
        vi, vj = p[i]/mass, p[j]/mass
        rrel = posi-posj  # vecteur pour la distance entre les centres des 2 sphères
        vrel = vj-vi   # vecteur pour la différence de vitesse entre les 2 sphères

        # exclusion de cas où il n'y a pas de changements à faire
        if np.sum((vrel)**2) == 0: continue  # exactly same velocities si et seulement si le vecteur vrel devient nul, la trajectoire des 2 sphères continue alors côte à côte
        if np.sum((pos[i] - pos[j])**2) > Ratom**2: continue  # one atom went all the way through another, la collision a été "manquée" à l'intérieur du pas deltax

        # calcule la distance et temps d'interpénétration des sphères dures qui ne doit pas se produire dans ce modèle
        dx = np.dot(rrel, vrel)/np.sqrt(np.sum(vrel**2))        # length of the component of rrel parallel to vrel
        dy = np.sum((np.cross(rrel, vrel)/np.sqrt(np.sum(vrel**2)))**2) # length of the component of rrel perpendicular to vrel
        alpha = np.asin(dy/(2*Ratom))       # alpha is the angle of the triangle composed of rrel, path of atom j, and a line from the center of atom i to the center of atom j where atome j hits atom i
        d = (2*Ratom)*np.cos(alpha)-dx      # distance traveled into the atom from first contact
        deltat = d/np.sqrt(np.sum(vrel**2)) # time spent moving from first contact to position inside atom

        #### CHANGE L'INTERPÉNÉTRATION DES SPHÈRES PAR LA CINÉTIQUE DE COLLISION ####
        posi, posj = posi-vi*deltat, posj-vj*deltat # back up to contact configuration
        pcmi, pcmj = p[i]-mass*vcm, p[j]-mass*vcm   # transform momenta to center-of-momentum frame
        pcmi, pcmj = pcmi-2*rrel*np.sum(pcmi*rrel)/np.sum(rrel**2), pcmj-2*rrel*np.sum(pcmj*rrel)/np.sum(rrel**2) # bounce in center-of-momentum frame
        p[i], p[j] = pcmi+mass*vcm, pcmj+mass*vcm # transform momenta back to lab frame
        pos[i], pos[j] = posi+(p[i]/mass)*deltat, posj+(p[j]/mass)*deltat # move forward deltat in time, ramenant au même temps où sont rendues les autres sphères dans l'itération

targeted_particle = follow_particle(targeted_particle, target_particle_index, step, dt, clean=True)
