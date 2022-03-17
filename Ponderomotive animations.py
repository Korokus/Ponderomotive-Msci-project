from enum import auto
import math
import numpy as np
import random as rand
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import os
from matplotlib import animation
from mpmath import *
#How the enhancement varies with different resonance conditions

#Add elipsoid
#Add damping -have to figure it out.
#Write a project plan
#Consider only ellipsoids of revolution

#Input parameters in SI units
#sphere radius 50nm
sphere_radius = 50e-9
#wavelength 500nm
wavelength = 500e-9
#Amplitude V/m
#1mw comfortable level -> 1 micron assuming some power density -> known flux
amplitude = 1e5
sphere_permittivity = -2.001
surrounding_permittivity = 1
c = 299792458
e = 1.60217662e-19
e_0 = 8.85418782e-12
direction = [1, 1, 0]
timesteps = 50000
electron_mass = 9.10938356e-31
omega = 2 * math.pi * c / wavelength
#setting phi and theta and r.  theta [0,pi] phi [0,2pi]
phi = 0
theta = math.pi/4

#0.1 nm away from the sphere radius
r = sphere_radius + 0.0000000001
phi = phi % math.pi
theta = theta % (2*math.pi)
trajectories = 20

#setting the radius of the radius for the ellisoid axis here we only consider where R_1 = R_2 or R_2 = R_3
#define R_1 and R_2 to be in the x-y plane and R_3 corresponds to the z plane
#only define R_1 >= R_2 >= R_3
R_1 =  10
R_2 =  1
R_3 =  1


def main():
    timestep = interpolate()
    time = 0
    incident_field = np.zeros((timesteps, 3))
    #stores the incident field in an array
    for t in range(timesteps):
        Ex = amplitude * math.cos(omega*time) * direction[0]
        Ey = amplitude * math.cos(omega*time) * direction[1]
        Ez = amplitude * math.cos(omega*time) * direction[2]
        incident_field[t][0] = Ex
        incident_field[t][1] = Ey
        incident_field[t][2] = Ez
        time += timestep
    #generating coordinate points for electron
    electron_position = np.zeros((timesteps, 3))
    electron_velocity = np.zeros((timesteps, 3))
    x = r * math.cos(theta)
    y = r * math.cos(phi) * math.sin(theta)
    z = r * math.sin(phi) * math.sin(theta)
    electron_position[0][0] = x
    electron_position[0][1] = y
    electron_position[0][2] = z
    graph_field = []
    #calculating the electric field at the electron position and updating the position of the electron
    for t in range(timesteps - 1):
        current_field = incident_field[t]
        field = electric_field(
            electron_position[t][0], electron_position[t][1], electron_position[t][2], current_field)
        graph_field.append(field[0])
        if t == 0:
            electron_position[1][0] = electron_position[0][0] + \
                (field[0] * e/electron_mass) * (timestep**2)
            electron_position[1][1] = electron_position[0][1] + \
                (field[1] * e/electron_mass) * (timestep**2)
            electron_position[1][2] = electron_position[0][2] + \
                (field[2] * e/electron_mass) * (timestep**2)
        else:
            electron_position[t+1][0] = (2 * electron_position[t][0]) - \
                                         electron_position[t-1][0] + (timestep**2) * (
                                             field[0] * e/electron_mass)
            electron_position[t+1][1] = (2 * electron_position[t][1]) - \
                electron_position[t-1][1] + \
                (timestep**2) * (field[1] * e/electron_mass)
            electron_position[t+1][2] = (2 * electron_position[t][2]) - \
                electron_position[t-1][2] + \
                (timestep**2) * (field[2] * e/electron_mass)
    #calculating the velocity using the verlet algorithm
    for t in range(timesteps - 1):
        if t == 0:
            electron_velocity[0][0] = 0
            electron_velocity[0][1] = 0
            electron_velocity[0][2] = 0
        else:
            electron_velocity[t][0] = (
                electron_position[t+1][0] - electron_position[t-1][0]) / (2 * timestep)
            electron_velocity[t][1] = (
                electron_position[t+1][1] - electron_position[t-1][1]) / (2 * timestep)
            electron_velocity[t][2] = (
                electron_position[t+1][2] - electron_position[t-1][2]) / (2 * timestep)

    timestep_graph(electron_position, timestep, graph_field, electron_velocity)
    #Creating a grid of x,y plots
    E2 = np.vectorize(E)
    nx, ny = 100, 100  # level of discretisation
    xr = np.linspace(-5*sphere_radius, 5*sphere_radius, nx)
    yr = np.linspace(-5*sphere_radius, 5*sphere_radius, ny)
    X, Y = np.meshgrid(xr, yr)
    plot_Ex = np.zeros((nx, ny))
    plot_Ey = np.zeros((nx, ny))
    f1 = incident_field[t][0]
    ex, ey = E2(x=X, y=Y, f1=f1)
    plot_Ex = ex
    plot_Ey = ey
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_xlim(-5*sphere_radius, 5*sphere_radius)
    ax.set_ylim(-5*sphere_radius, 5*sphere_radius)
    ax.set_aspect('equal')
    ax.add_artist(Circle((0, 0), sphere_radius, fill=False))
    Q = ax.quiver(xr, yr, plot_Ex, plot_Ey, pivot='mid')
    anim = animation.FuncAnimation(fig, update_quiver, fargs=(Q, X, Y, incident_field, E2),
                                   interval=1000, blit=False)
    fig.tight_layout()
    # saving to mp4 using ffmpeg writer
    anim.save('Electric field.gif', fps=60, dpi=300)
    plt.close()
    E2 = np.vectorize(E)
    nx, ny = 1000, 1000  # level of discretisation
    xr = np.linspace(-5*sphere_radius, 5*sphere_radius, nx)
    yr = np.linspace(-5*sphere_radius, 5*sphere_radius, ny)
    X, Y = np.meshgrid(xr, yr)
    plot_Ex = np.zeros((nx, ny))
    plot_Ey = np.zeros((nx, ny))
    f1 = incident_field[t][0]
    ex, ey = E2(x=X, y=Y, f1=f1)
    plot_Ex = ex
    plot_Ey = ey
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_xlim(-5*sphere_radius,5*sphere_radius)
    ax.set_ylim(-5*sphere_radius,5*sphere_radius)
    ax.set_aspect('equal')
    # Plot the streamlines with an appropriate colormap and arrow style
    color = 2 * np.log(np.hypot(plot_Ex, plot_Ey))
    ax.streamplot(xr, yr, plot_Ex, plot_Ey, color=color, linewidth=0.2, cmap=plt.cm.inferno,
    density=5, arrowstyle='->', arrowsize=0.5)
    ax.add_artist(Circle((0, 0), sphere_radius, fill=False))
    plt.savefig("electric field lines streamplot",dpi=300, bbox_inches='tight')
    nx, ny = 1000, 1000  # level of discretisation
    xr = np.linspace(-5*sphere_radius, 5*sphere_radius, nx)
    yr = np.linspace(-5*sphere_radius, 5*sphere_radius, ny)
    X, Y = np.meshgrid(xr, yr)
    fig = plt.figure()
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.xlim(-3*sphere_radius, 3*sphere_radius)
    plt.ylim(-3*sphere_radius, 3*sphere_radius)
    enhancement2 = np.vectorize(enhancement)
    Z = enhancement2(x=X, y=Y)
    plt.pcolormesh(xr, yr, Z)
    plt.colorbar()
    angles = np.linspace(-math.pi/2,math.pi/2, trajectories)
    for i in range(trajectories -1):
        phi_trajectory = 0
        theta_trajectory = angles[i]
        r_trajectory =  sphere_radius + 0.0000000001
        phi_trajectory = phi_trajectory % math.pi
        theta_trajectory = theta_trajectory % (2*math.pi)
        electron_positions = np.zeros((timesteps,3))
        graph_position = []
        electron_velocity_position = np.zeros((timesteps,3))
        x = r_trajectory * math.cos(theta_trajectory)
        y = r_trajectory * math.cos(phi_trajectory) * math.sin(theta_trajectory)
        z = r_trajectory * math.sin(phi_trajectory) * math.sin(theta_trajectory)
        electron_positions[0][0] = x
        electron_positions[0][1] = y
        electron_positions[0][2] = z
        xpos = []
        ypos = []
        xvel = []
        yvel = []
        electron_positions, graph_position, electron_velocity_position = positions(electron_positions, graph_position, incident_field, electron_velocity_position, timestep)
        for t in range(timesteps-1):
            xpos.append(electron_positions[t][0])
            ypos.append(electron_positions[t][1])
            xvel.append(electron_velocity_position[t][0])
            yvel.append(electron_velocity_position[t][1])
        plt.plot(xpos, ypos)
    plt.savefig("electric field lines density colormesh",
                dpi=300, bbox_inches='tight')
    plt.close()
    fig = plt.figure()
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.xlim(-3*sphere_radius, 3*sphere_radius)
    plt.ylim(-3*sphere_radius, 3*sphere_radius)
    enhancement2 = np.vectorize(enhancement_X)
    Z = enhancement2(x=X, y=Y)
    plt.pcolormesh(xr, yr, Z)
    plt.colorbar()
    plt.savefig("electric field lines density colormesh X direction",
                dpi=300, bbox_inches='tight')
    plt.close()
    fig = plt.figure()
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.xlim(-3*sphere_radius, 3*sphere_radius)
    plt.ylim(-3*sphere_radius, 3*sphere_radius)
    enhancement2 = np.vectorize(enhancement_Y)
    Z = enhancement2(x=X, y=Y)
    plt.pcolormesh(xr, yr, Z)
    plt.colorbar()
    plt.savefig("electric field lines density colormesh Y direction",
                dpi=300, bbox_inches='tight')
    plt.close()


    #Code for ellipsoidal particles

    #Finding the L_factors for the different axis. 
    L_factor_R_1,L_factor_R_2,L_factor_R_3 = L_factor(R_1,R_2,R_3)
    ellipsoid_field2 = np.vectorize(E_ellipsoid_outside)
    timesteps2 = 20000
    electron_position_ellipse = np.zeros((timesteps2, 3))
    electron_velocity_ellipse = np.zeros((timesteps2, 3))
    x = r * math.cos(theta)
    y = r * math.cos(phi) * math.sin(theta)
    z = r * math.sin(phi) * math.sin(theta)
    electron_position_ellipse[0][0] = x
    electron_position_ellipse[0][1] = y
    electron_position_ellipse[0][2] = z
    for t in range(timesteps2 - 1):
        print(t)
        current_field = incident_field[t]
        field = E_ellipsoid_outside(
            electron_position_ellipse[t][0], electron_position[t][1], electron_position[t][2], L_factor_R_1, L_factor_R_2, L_factor_R_3, current_field)
        graph_field.append(field[0])
        if t == 0:
            electron_position_ellipse[1][0] = electron_position_ellipse[0][0] + \
                (field[0] * e/electron_mass) * (timestep**2)
            electron_position_ellipse[1][1] = electron_position_ellipse[0][1] + \
                (field[1] * e/electron_mass) * (timestep**2)
            electron_position_ellipse[1][2] = electron_position_ellipse[0][2] + \
                (field[2] * e/electron_mass) * (timestep**2)
        else:
            electron_position_ellipse[t+1][0] = (2 * electron_position_ellipse[t][0]) - \
                                         electron_position_ellipse[t-1][0] + (timestep**2) * (
                                             field[0] * e/electron_mass)
            electron_position_ellipse[t+1][1] = (2 * electron_position_ellipse[t][1]) - \
                electron_position_ellipse[t-1][1] + \
                (timestep**2) * (field[1] * e/electron_mass)
            electron_position_ellipse[t+1][2] = (2 * electron_position_ellipse[t][2]) - \
                electron_position_ellipse[t-1][2] + \
                (timestep**2) * (field[2] * e/electron_mass)
    
    #calculating the velocity using the verlet algorithm
    for t in range(timesteps2 - 1):
        if t == 0:
            electron_velocity_ellipse[0][0] = 0
            electron_velocity_ellipse[0][1] = 0
            electron_velocity_ellipse[0][2] = 0
        else:
            electron_velocity[t][0] = (
                electron_position_ellipse[t+1][0] - electron_position_ellipse[t-1][0]) / (2 * timestep)
            electron_velocity_ellipse[t][1] = (
                electron_position_ellipse[t+1][1] - electron_position_ellipse[t-1][1]) / (2 * timestep)
            electron_velocity_ellipse[t][2] = (
                electron_position_ellipse[t+1][2] - electron_position_ellipse[t-1][2]) / (2 * timestep)
    timestep_graph_ellipse(electron_position_ellipse,timesteps2,electron_velocity_ellipse,timestep)



def inside_check(x, y, z):
    r = math.sqrt((x**2) + (y**2) + (z**2))
    if r < sphere_radius:
        return True
    else:
        return False
#Returns the electric field for the quiver plot


def E(x, y, f1):
    z = 0
    r = math.sqrt((x**2) + (y**2))

    if inside_check(x, y, z) == False:
        Ex = (1 + ((sphere_permittivity - surrounding_permittivity)/(sphere_permittivity + 2
              * surrounding_permittivity)) * ((sphere_radius**3)/r**5) * (2 * (x**2) - (y**2) - (z**2))) * f1
        Ey = ((sphere_permittivity - surrounding_permittivity)/(sphere_permittivity
              + 2 * surrounding_permittivity)) * ((sphere_radius**3)/r**5) * (3*x*y) * f1
        return Ex, Ey
    if inside_check(x, y, z) == True:
        Ex = f1 * (3*surrounding_permittivity
                   / (sphere_permittivity + 2 * surrounding_permittivity))
        Ey = 0
        return Ex, Ey

#Returns the electric field for the electron


def electric_field(x, y, z, current_field):
    r = math.sqrt((x**2) + (y**2) + (z**2))
    if inside_check(x, y, z) == False:
        Ex = (1 + ((sphere_permittivity - surrounding_permittivity)/(sphere_permittivity + 2 * surrounding_permittivity))
              * ((sphere_radius**3)/r**5) * (2 * (x**2) - (y**2) - (z**2))) * current_field[0]
        Ey = ((sphere_permittivity - surrounding_permittivity)/(sphere_permittivity + 2
              * surrounding_permittivity)) * ((sphere_radius**3)/r**5) * (3*x*y) * current_field[0]
        Ez = ((sphere_permittivity - surrounding_permittivity)/(sphere_permittivity + 2
              * surrounding_permittivity)) * ((sphere_radius**3)/r**5) * (3*x*z) * current_field[0]
        return[Ex, Ey, Ez]
    if inside_check(x, y, z) == True:
        Ex = current_field[0] * (3*surrounding_permittivity
                                 / (sphere_permittivity + 2 * surrounding_permittivity))
        Ey = 0
        Ez = 0
        return[Ex, Ey, Ez]
#function that plots different graphs


def timestep_graph(electron_position, timestep, graph_field, electron_velocity):
    timestep_for_graph = []
    xpos = []
    ypos = []
    zpos = []
    xvel = []
    yvel = []
    zvel = []
    temp = 0
    for t in range(timesteps-1):
        timestep_for_graph.append(temp)
        temp += timestep
        xpos.append(electron_position[t][0])
        ypos.append(electron_position[t][1])
        zpos.append(electron_position[t][2])
        xvel.append(electron_velocity[t][0])
        yvel.append(electron_velocity[t][1])
        zvel.append(electron_velocity[t][2])
    plt.figure(0)
    plt.plot(timestep_for_graph, xpos)
    plt.savefig("Xposition of electron.png", dpi=300, bbox_inches='tight')
    plt.close()
    plt.figure(1)
    plt.plot(timestep_for_graph, ypos)
    plt.savefig("Yposition of electron.png", dpi=300, bbox_inches='tight')
    plt.close()
    plt.figure(2)
    plt.plot(timestep_for_graph, zpos)
    plt.savefig("Zposition of electron.png", dpi=300, bbox_inches='tight')
    plt.close()
    plt.figure(3)
    plt.plot(timestep_for_graph, graph_field)
    plt.savefig("electric field amplitude.png", dpi=300, bbox_inches='tight')
    plt.close()
    plt.plot(timestep_for_graph, xvel)
    plt.savefig("Xvelocity of electron.png", dpi=300, bbox_inches='tight')
    plt.close()
    plt.plot(timestep_for_graph, yvel)
    plt.savefig("Yvelocity of electron.png", dpi=300, bbox_inches='tight')
    plt.close()
    plt.plot(timestep_for_graph, zvel)
    plt.savefig("Zvelocity of electron.png", dpi=300, bbox_inches='tight')
    plt.close()
    plt.plot(xpos, ypos)
    plt.savefig("X-Y trajectory of electron.png", dpi=300, bbox_inches='tight')
    plt.close()
    plt.plot(xpos, zpos)
    plt.savefig("X-Z trajectory of electron.png", dpi=300, bbox_inches='tight')
    plt.close()

def timestep_graph_ellipse(electron_position, timesteps2, electron_velocity,timestep):
    timestep_for_graph = []
    xpos = []
    ypos = []
    zpos = []
    xvel = []
    yvel = []
    zvel = []
    temp = 0
    for t in range(timesteps2-1):
        timestep_for_graph.append(temp)
        temp += timestep
        xpos.append(electron_position[t][0])
        ypos.append(electron_position[t][1])
        zpos.append(electron_position[t][2])
        xvel.append(electron_velocity[t][0])
        yvel.append(electron_velocity[t][1])
        zvel.append(electron_velocity[t][2])
    plt.plot(timestep_for_graph, xpos)
    plt.savefig("Xposition of electron ellipse.png", dpi=300, bbox_inches='tight')
    plt.close()
    plt.plot(timestep_for_graph, ypos)
    plt.savefig("Yposition of electron ellipse.png", dpi=300, bbox_inches='tight')
    plt.close()
    plt.plot(timestep_for_graph, zpos)
    plt.savefig("Zposition of electron ellipse.png", dpi=300, bbox_inches='tight')
    plt.close()

#function calculates the timestep
def interpolate():
    hertz = c/wavelength
    timestep = 1/(hertz * 10)
    return timestep

#quiverplot animation


def update_quiver(num, Q, X, Y, incident_field, E2):
    f1 = incident_field[num][0]
    ex, ey = E2(x=X, y=Y, f1=f1)
    Q.set_UVC(ex, ey)
    return Q
#calculates the enhancement of the field


def enhancement(x, y):
    z = 0
    r = math.sqrt((x**2) + (y**2))
    if inside_check(x, y, z) == False:
        Ex = (1 + ((sphere_permittivity - surrounding_permittivity)/(sphere_permittivity + 2
              * surrounding_permittivity)) * ((sphere_radius**3)/r**5) * (2 * (x**2) - (y**2) - (z**2)))
        Ey = ((sphere_permittivity - surrounding_permittivity)/(sphere_permittivity
              + 2 * surrounding_permittivity)) * ((sphere_radius**3)/r**5) * (3*x*y)
        return math.sqrt((Ex**2) + (Ey**2))
    if inside_check(x, y, z) == True:
        Ex = (3*surrounding_permittivity
              / (sphere_permittivity + 2 * surrounding_permittivity))
        Ey = 0
        return math.sqrt((Ex**2) + (Ey**2))


def enhancement_X(x, y):
    z = 0
    r = math.sqrt((x**2) + (y**2))
    if inside_check(x, y, z) == False:
        Ex = (1 + ((sphere_permittivity - surrounding_permittivity)/(sphere_permittivity + 2
              * surrounding_permittivity)) * ((sphere_radius**3)/r**5) * (2 * (x**2) - (y**2) - (z**2)))
        Ey = 0
        return math.sqrt((Ex**2) + (Ey**2))
    if inside_check(x, y, z) == True:
        Ex = (3*surrounding_permittivity
              / (sphere_permittivity + 2 * surrounding_permittivity))
        Ey = 0
        return math.sqrt((Ex**2) + (Ey**2))


def enhancement_Y(x, y):
    z = 0
    r = math.sqrt((x**2) + (y**2))
    if inside_check(x, y, z) == False:
        Ex = 0
        Ey = ((sphere_permittivity - surrounding_permittivity)/(sphere_permittivity
              + 2 * surrounding_permittivity)) * ((sphere_radius**3)/r**5) * (3*x*y)
        return math.sqrt((Ex**2) + (Ey**2))
    if inside_check(x, y, z) == True:
        Ex = 0
        Ey = 0
        return math.sqrt((Ex**2) + (Ey**2))

def positions(electron_position, graph_field, incident_field, electron_velocity, timestep):
    for t in range(timesteps - 1):
        current_field = incident_field[t]
        field = electric_field(
            electron_position[t][0], electron_position[t][1], electron_position[t][2], current_field)
        graph_field.append(field[0])
        if t == 0:
            electron_position[1][0] = electron_position[0][0] + \
                (field[0] * e/electron_mass) * (timestep**2)
            electron_position[1][1] = electron_position[0][1] + \
                (field[1] * e/electron_mass) * (timestep**2)
            electron_position[1][2] = electron_position[0][2] + \
                (field[2] * e/electron_mass) * (timestep**2)
        else:
            electron_position[t+1][0] = (2 * electron_position[t][0]) - \
                                         electron_position[t-1][0] + (timestep**2) * (
                                             field[0] * e/electron_mass)
            electron_position[t+1][1] = (2 * electron_position[t][1]) - \
                electron_position[t-1][1] + \
                (timestep**2) * (field[1] * e/electron_mass)
            electron_position[t+1][2] = (2 * electron_position[t][2]) - \
                electron_position[t-1][2] + \
                (timestep**2) * (field[2] * e/electron_mass)
    

    #calculating the velocity using the verlet algorithm
    for t in range(timesteps - 1):
        if t == 0:
            electron_velocity[0][0] = 0
            electron_velocity[0][1] = 0
            electron_velocity[0][2] = 0
        else:
            electron_velocity[t][0] = (
                electron_position[t+1][0] - electron_position[t-1][0]) / (2 * timestep)
            electron_velocity[t][1] = (
                electron_position[t+1][1] - electron_position[t-1][1]) / (2 * timestep)
            electron_velocity[t][2] = (
                electron_position[t+1][2] - electron_position[t-1][2]) / (2 * timestep)
    return electron_position, graph_field, electron_velocity

def L_factor(R_1,R_2,R_3):
    if R_1 == R_2 > R_3:
        L_factor_R_3 = prolate_definite(R_1,R_3)
        L_factor_R_1 = 0.5 * (1-L_factor_R_3)
        L_factor_R_2 = 0.5 * (1-L_factor_R_3)
        
    if R_1 > R_2 == R_3:
        L_factor_R_1 = oblate_definite(R_1,R_3)
        L_factor_R_2 = 0.5 * (1 - L_factor_R_1)
        L_factor_R_3 = 0.5 * (1 - L_factor_R_1)

    return L_factor_R_1,L_factor_R_2,L_factor_R_3

        

#Returns the electric field for the ellispoidal particle 
def E_ellipsoid_inside(L_factor_R_1, L_factor_R_2,L_factor_R_3,incident_field):
    ellipsoid_field = np.zeros((timesteps, 3))
    for t in range(timesteps -1):
        current_field = incident_field[t]
        P_x = e_0 * (sphere_permittivity - surrounding_permittivity) * ((surrounding_permittivity) * current_field[0])/(surrounding_permittivity + L_factor_R_1 (sphere_permittivity - surrounding_permittivity))
        P_y = e_0 * (sphere_permittivity - surrounding_permittivity) * ((surrounding_permittivity) * current_field[1])/(surrounding_permittivity + L_factor_R_2 (sphere_permittivity - surrounding_permittivity))
        P_z = e_0 * (sphere_permittivity - surrounding_permittivity) * ((surrounding_permittivity) * current_field[2])/(surrounding_permittivity + L_factor_R_3 (sphere_permittivity - surrounding_permittivity))
        ellipsoid_field[t][0] = current_field[0] - L_factor_R_1* P_x/(e_0 * surrounding_permittivity)
        ellipsoid_field[t][1] = current_field[1] - L_factor_R_2* P_y/(e_0 * surrounding_permittivity)
        ellipsoid_field[t][2] = current_field[2] - L_factor_R_3* P_z/(e_0 * surrounding_permittivity)
    return ellipsoid_field

def E_ellipsoid_outside(x,y,z,L_factor_R_1,L_factor_R_2,L_factor_R_3,current_field):
    prefactor = (R_1*R_2*R_3/2)* ((surrounding_permittivity - sphere_permittivity)/sphere_permittivity)
    xi = Xi(x,y,z)
    if R_1 == R_2 > R_3:
        integral_x = oblate_indefinite(xi,R_1,R_3)
        integral_y = oblate_indefinite(xi,R_1,R_3)
        integral_z = prolate_indefinite(xi,R_1,R_3)
    if R_1 > R_2 == R_3:
        integral_x = prolate_indefinite(xi,R_1,R_3)
        integral_y = prolate_indefinite(xi,R_1,R_3)
        integral_z = oblate_indefinite(xi,R_1,R_3)
    if ellipsoid_inside_check(x,y,z) == True:
        E_x = (1 + prefactor * integral_x)/(1 + prefactor*L_factor_R_1) * current_field[0]
        E_y = (1 + prefactor * integral_y)/(1 + prefactor*L_factor_R_2) * current_field[1]
        E_z = (1 + prefactor * integral_z)/(1 + prefactor*L_factor_R_3) * current_field[2]
    else:
        E_x = current_field[0] / (prefactor * L_factor_R_1)
        E_y = current_field[1] / (prefactor * L_factor_R_2)
        E_z = current_field[2] / (prefactor * L_factor_R_3)
    return [E_x,E_y,E_z]


def ellipsoid_inside_check(x,y,z):
    a = R_1 * 1e-9
    b = R_2 * 1e-9
    c = R_3 * 1e-9

    if (x**2)/(a**2) + (y**2)/(b**2) + (z**2)/(c**2) > 1:
        return True
    else:
        return False






def Xi(x,y,z):
    h = math.sqrt(R_1**2 - R_2**2)
    k = math.sqrt(R_1**2 - R_3**2)
    a_1 = -((x**2)+(y**2)+(z**2)+(h**2)+(k**2))
    a_2 = (x**2)*((h**2) + (k**2)) + (y**2)*(k**2) + (z**2)*(h**2) + (h**2)* (k**2)
    a_3 = (x**2)*(h**2)*(k**2)
    Q = ((a_1**2) - 3 * a_2) / 9
    R = ((9*a_1*a_2) - (27 * a_3) - (2 * a_1**3))/54
    theta = math.acos(R/math.sqrt(Q**3))
    xi = 2 * math.sqrt(Q) * math.cos (theta/3) - (a_1/3)
    return math.sqrt(xi)



def oblate_indefinite(xi,a,b):
    x = a**2
    y = b**2
    num_1 = math.sqrt(xi + y)
    num_2 = math.atan(math.sqrt(xi + y)/math.sqrt(x-y))
    denom_1 = (xi + x) * (xi - y)
    denom_2 = (x-y)**(3/2)
    return (num_1/denom_1) + (num_2/denom_2)

def oblate_definite(a,b):
    e = math.sqrt(1-b**2/a**2)
    prefactor = -1/((a**3)*e**3)
    num = 1+e 
    denom = 1-e 
    return prefactor *(2*e - math.log(num/denom))

def prolate_indefinite(xi,a,b):
    x = a**2
    y = b**2
    mp.dps = 25
    mp.pretty = True
    #note only defnined for abs (xi + y)/(y-x) < 1
    num_1 = 2 * hyp2f1(-0.5,1,0.5,(xi + y)/(y-x))
    denom_1 = math.sqrt(xi + y) * (x - y)
    return num_1/denom_1

def prolate_definite(a,c):
    x = a**2
    y = c**2
    num_1 = math.sqrt(x-y)/c
    num_2 = math.atan(math.sqrt(x-y)/c)
    denom_1 = (x-y)**(3/2)
    return 2 * (num_1 - num_2) / denom_1





if __name__ == '__main__':
    main()
