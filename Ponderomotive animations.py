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
direction = [1, 0, 0]
timesteps = 200000
electron_mass = 9.10938356e-31
omega = 2 * math.pi * c / wavelength
#setting phi and theta and r.  theta [0,pi] phi [0,2pi]
phi = 0
theta = -math.pi/3

#0.1 nm away from the sphere radius
r = sphere_radius + 0.0000000001
phi = phi % math.pi
theta = theta % (2*math.pi)
trajectories = 20

#setting the radius of the radius for the ellisoid axis here we only consider where R_1 = R_2 or R_2 = R_3
#define R_1 and R_2 to be in the x-y plane and R_3 corresponds to the z plane
#only define R_1 >= R_2 >= R_3
R_1 =  5
R_2 =  1
R_3 =  1

#setting the initial phase of the plane wave
phase = math.pi/6

def main():
    timestep = interpolate()
    time = 0
    incident_field = np.zeros((timesteps, 3))
    #stores the incident field in an array
    for t in range(timesteps):
        Ex = amplitude * math.cos(omega*time + phase) * direction[0]
        Ey = amplitude * math.cos(omega*time + phase) * direction[1]
        Ez = amplitude * math.cos(omega*time + phase) * direction[2]
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
    plt.close()
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
    angles = np.linspace(0,math.pi*2, trajectories)
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
        electron_positions, graph_position, electron_velocity_position = positions(electron_positions, graph_position, incident_field, electron_velocity_position, timestep,x,y,z,timesteps)
        for t in range(timesteps-1):
            xpos.append(electron_positions[t][0])
            ypos.append(electron_positions[t][1])
            xvel.append(electron_velocity_position[t][0])
            yvel.append(electron_velocity_position[t][1])
        plt.plot(xpos, ypos)
    plt.savefig("electric field lines density colormesh",
                dpi=300, bbox_inches='tight')
    plt.close()
    for i in range(trajectories -1):
        phi_trajectory = 0
        theta_trajectory = angles[i]
        r_trajectory =  sphere_radius + 0.0000000001
        phi_trajectory = phi_trajectory % math.pi
        theta_trajectory = theta_trajectory % (2*math.pi)
        electron_positions = np.zeros((20,3))
        graph_position = []
        electron_velocity_position = np.zeros((20,3))
        x = r_trajectory * math.cos(theta_trajectory)
        y = r_trajectory * math.cos(phi_trajectory) * math.sin(theta_trajectory)
        z = r_trajectory * math.sin(phi_trajectory) * math.sin(theta_trajectory)
        electron_positions[0][0] = x
        electron_positions[0][1] = y
        electron_positions[0][2] = z
        radius = math.sqrt(x**2 + y**2)
        dist = []
        time = []
        init_time = 0 
        electron_positions, graph_position, electron_velocity_position = positions(electron_positions, graph_position, incident_field, electron_velocity_position, timestep,x,y,z,20)
        for t in range(20):
            time.append(init_time)
            dist.append(math.sqrt(electron_positions[t][0]**2 + electron_positions[t][1]**2) - radius)
            init_time += timestep
        plt.plot(time,dist)
    plt.savefig("electron sphere trajectories",
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
    ellipsoid_permittivity = (((-1/L_factor_R_1) + 1) * surrounding_permittivity) + 0.0001
    free_space = e_0
    electron_position_ellipse = np.zeros((timesteps, 3))
    electron_velocity_ellipse = np.zeros((timesteps, 3))
    a = R_1 * 1e-8
    b = R_2 * 1e-8
    if 0 < theta < math.pi/4   or math.pi*3/4 < theta < 2*math.pi:
        x = (a*b) / (math.sqrt(b**2) + (a**2) * (math.tan(theta)**2))
        y = (a*b)*(math.tan(theta)) / (math.sqrt(b**2) + (a**2) * (math.tan(theta)**2))
    else: 
        x = -(a*b) / (math.sqrt(b**2) + (a**2) * (math.tan(theta)**2))
        y = -(a*b)*(math.tan(theta)) / (math.sqrt(b**2) + (a**2) * (math.tan(theta)**2))
    #z = r * math.sin(phi) * math.sin(theta)
    z = 0
    electron_position_ellipse[0][0] = x
    electron_position_ellipse[0][1] = y
    electron_position_ellipse[0][2] = z
    for t in range(timesteps - 1):
        current_field = incident_field[t]
        field = E_ellipsoid_outside(
            electron_position_ellipse[t][0], electron_position_ellipse[t][1], electron_position_ellipse[t][2], L_factor_R_1,current_field, ellipsoid_permittivity,free_space)
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
    for t in range(timesteps - 1):
        if t == 0:
            electron_velocity_ellipse[0][0] = 0
            electron_velocity_ellipse[0][1] = 0
            electron_velocity_ellipse[0][2] = 0
        else:
            electron_velocity_ellipse[t][0] = (
                electron_position_ellipse[t+1][0] - electron_position_ellipse[t-1][0]) / (2 * timestep)
            electron_velocity_ellipse[t][1] = (
                electron_position_ellipse[t+1][1] - electron_position_ellipse[t-1][1]) / (2 * timestep)
            electron_velocity_ellipse[t][2] = (
                electron_position_ellipse[t+1][2] - electron_position_ellipse[t-1][2]) / (2 * timestep)
    timestep_graph_ellipse(electron_position_ellipse,timesteps,electron_velocity_ellipse,timestep)
    E_ellipsoid2 = np.vectorize(E_ellipsoid)
    nx, ny = 1000, 1000  # level of discretisation
    xr = np.linspace(-5*sphere_radius, 5*sphere_radius, nx)
    yr = np.linspace(-5*sphere_radius, 5*sphere_radius, ny)
    X, Y = np.meshgrid(xr, yr)
    plot_Ex = np.zeros((nx, ny))
    plot_Ey = np.zeros((nx, ny))
    f1 = incident_field[0][0]
    f2 = incident_field[0][1]
    f3 = incident_field[0][2]
    ex, ey = E_ellipsoid2(x=X, y=Y, f1 = f1, L_factor_R_1=L_factor_R_1,ellipsoid_permittivity = ellipsoid_permittivity)
    plot_Ex = ex
    plot_Ey = ey
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_xlim(-5*sphere_radius,5*sphere_radius)
    ax.set_ylim(-5*sphere_radius,5*sphere_radius)
    ax.set_aspect('equal')
    color = 2 * np.log(np.hypot(abs(plot_Ex), abs(plot_Ey)))
    ax.streamplot(xr, yr, plot_Ex, plot_Ey, color=color, linewidth=0.2, cmap=plt.cm.inferno,
    density=5, arrowstyle='->', arrowsize=0.5)
    plt.savefig("electric field lines streamplot ellipsoid",dpi=300, bbox_inches='tight')
    plt.close()
    nx, ny = 1000, 1000  # level of discretisation
    xr = np.linspace(-5*sphere_radius, 5*sphere_radius, nx)
    yr = np.linspace(-5*sphere_radius, 5*sphere_radius, ny)
    X, Y = np.meshgrid(xr, yr)
    fig = plt.figure()
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.xlim(-3*sphere_radius, 3*sphere_radius)
    plt.ylim(-3*sphere_radius, 3*sphere_radius)
    enhancement_ellipsoid2 = np.vectorize(enhancement_ellipsoid)
    Z = enhancement_ellipsoid2(x=X, y=Y,f1 = f1, L_factor_R_1 = L_factor_R_1, ellipsoid_permittivity = ellipsoid_permittivity)
    plt.pcolormesh(xr, yr, Z)
    plt.colorbar()
    angles = np.linspace(-math.pi/2,math.pi/2, trajectories)
    for i in range(trajectories -1):
        phi_trajectory = 0
        theta_trajectory = angles[i]
        phi_trajectory = phi_trajectory % math.pi
        theta_trajectory = theta_trajectory % (2*math.pi)
        electron_positions = np.zeros((timesteps,3))
        graph_position = []
        electron_velocity_position = np.zeros((timesteps,3))
        if 0 < theta < math.pi/4   or math.pi*3/4 < theta < 2*math.pi:
            x = (a*b) / (math.sqrt(b**2) + (a**2) * (math.tan(theta)**2))
            y = (a*b)*(math.tan(theta)) / (math.sqrt(b**2) + (a**2) * (math.tan(theta)**2))
        else: 
            x = -(a*b) / (math.sqrt(b**2) + (a**2) * (math.tan(theta)**2))
            y = -(a*b)*(math.tan(theta)) / (math.sqrt(b**2) + (a**2) * (math.tan(theta)**2))
        z = 0
        electron_positions[0][0] = x
        electron_positions[0][1] = y
        electron_positions[0][2] = z
        xpos = []
        ypos = []
        xvel = []
        yvel = []
        electron_positions, graph_position, electron_velocity_position = positions_ellipsoid(electron_positions, graph_position, incident_field, electron_velocity_position, timestep,L_factor_R_1,ellipsoid_permittivity,free_space)
        for t in range(timesteps-1):
            xpos.append(electron_positions[t][0])
            ypos.append(electron_positions[t][1])
            xvel.append(electron_velocity_position[t][0])
            yvel.append(electron_velocity_position[t][1])
        plt.plot(xpos, ypos)
    plt.savefig("electric field lines density colormesh ellipsoid",
                dpi=300, bbox_inches='tight')
    print("done")

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
        #print(Ex,Ey)
        return Ex, Ey
    if inside_check(x, y, z) == True:
        Ex = f1 * (3*surrounding_permittivity
                   / (sphere_permittivity + 2 * surrounding_permittivity))
        Ey = 0
        #print(Ex,Ey)
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

def timestep_graph_ellipse(electron_position, timesteps, electron_velocity,timestep):
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
    plt.plot(timestep_for_graph, xpos)
    plt.savefig("Xposition of electron ellipse.png", dpi=300, bbox_inches='tight')
    plt.close()
    plt.plot(timestep_for_graph, ypos)
    plt.savefig("Yposition of electron ellipse.png", dpi=300, bbox_inches='tight')
    plt.close()
    plt.plot(timestep_for_graph, zpos)
    plt.savefig("Zposition of electron ellipse.png", dpi=300, bbox_inches='tight')
    plt.close()
    plt.plot(xpos, ypos)
    plt.savefig("X-Y trajectory of electron ellipse.png", dpi=300, bbox_inches='tight')
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

def positions(electron_position, graph_field, incident_field, electron_velocity, timestep,x,y,z,cycles):
    for t in range(cycles - 1):
        current_field = incident_field[t]
        field = electric_field(
            electron_position[t][0], electron_position[t][1], electron_position[t][2], current_field)
        if inside_check(electron_position[t][0],electron_position[t][1],electron_position[t][2]) == True:
            electron_position[t+1][0] = x * 0.95
            electron_position[t+1][1] = y * 0.95
            electron_position[t+1][2] = z * 0.95
        elif t == 0:
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
    for t in range(cycles - 1):
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

def positions_ellipsoid(electron_position_ellipse, graph_field, incident_field, electron_velocity_ellipse, timestep,L_factor_R_1,ellipsoid_permittivity,free_space):
    for t in range(timesteps - 1):
        current_field = incident_field[t]
        field = E_ellipsoid_outside(
            electron_position_ellipse[t][0], electron_position_ellipse[t][1], electron_position_ellipse[t][2], L_factor_R_1,current_field, ellipsoid_permittivity,free_space)
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
    for t in range(timesteps - 1):
        if t == 0:
            electron_velocity_ellipse[0][0] = 0
            electron_velocity_ellipse[0][1] = 0
            electron_velocity_ellipse[0][2] = 0
        else:
            electron_velocity_ellipse[t][0] = (
                electron_position_ellipse[t+1][0] - electron_position_ellipse[t-1][0]) / (2 * timestep)
            electron_velocity_ellipse[t][1] = (
                electron_position_ellipse[t+1][1] - electron_position_ellipse[t-1][1]) / (2 * timestep)
            electron_velocity_ellipse[t][2] = (
                electron_position_ellipse[t+1][2] - electron_position_ellipse[t-1][2]) / (2 * timestep)
    return electron_position_ellipse, graph_field, electron_velocity_ellipse


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

def E_ellipsoid(x,y,f1,L_factor_R_1,ellipsoid_permittivity):
    prefactor = (R_1*R_2*R_3/2)* ((surrounding_permittivity - ellipsoid_permittivity)/ellipsoid_permittivity)
    V = (R_1*R_2*R_3) * 1e-24
    z = 0
    r = math.sqrt(x**2 + y**2 + z**2)
    P_x = f1 / (L_factor_R_1 + 1/(ellipsoid_permittivity + surrounding_permittivity))
    if ellipsoid_inside_check(x,y,z) == False:
        E_x = -(f1 - L_factor_R_1* P_x)
        E_y = 0
    else:
        E_x = f1 * (1 + surrounding_permittivity*(ellipsoid_permittivity - surrounding_permittivity) * ((surrounding_permittivity) )/(surrounding_permittivity + L_factor_R_1 * (ellipsoid_permittivity - surrounding_permittivity)) * (V/r**5) * (2 * (x**2) - (y**2) - (z**2)))
        E_y = f1 * (surrounding_permittivity * (ellipsoid_permittivity - surrounding_permittivity) * ((surrounding_permittivity) )/(surrounding_permittivity + L_factor_R_1 * (ellipsoid_permittivity - surrounding_permittivity)) * (V/r**5) * (3*x*y))
    return E_x,E_y
    

def E_ellipsoid_outside(x,y,z,L_factor_R_1,current_field,ellipsoid_permittivity,free_space):
    #prefactor = (R_1*R_2*R_3/2)* ((surrounding_permittivity - ellipsoid_permittivity)/ellipsoid_permittivity)
    V = (R_1*R_2*R_3)*1e-24
    r = math.sqrt(x**2 + y**2 + z**2)
    P_x = free_space * (sphere_permittivity - surrounding_permittivity) * ((surrounding_permittivity) * current_field[0])/(surrounding_permittivity + L_factor_R_1 * (sphere_permittivity - surrounding_permittivity))
    if ellipsoid_inside_check(x,y,z) == False:
        E_x = (current_field[0] - L_factor_R_1* P_x/(e_0 * surrounding_permittivity))
        E_y = 0
        E_z = 0
    else:
        E_x = current_field[0] * (1 + (ellipsoid_permittivity - surrounding_permittivity) * ((surrounding_permittivity))/(surrounding_permittivity + L_factor_R_1 * (ellipsoid_permittivity - surrounding_permittivity)) * (V/r**5) * (2 * (x**2) - (y**2) - (z**2)))
        E_y = current_field[0] * (ellipsoid_permittivity - surrounding_permittivity) * ((surrounding_permittivity) )/(surrounding_permittivity + L_factor_R_1 * (ellipsoid_permittivity - surrounding_permittivity)) * (V/r**5) * (3*x*y)
        E_z = current_field[0] * (ellipsoid_permittivity - surrounding_permittivity) * ((surrounding_permittivity) )/(surrounding_permittivity + L_factor_R_1 * (ellipsoid_permittivity - surrounding_permittivity)) * (V/r**5) * (3*z*y)
    return [E_x,E_y,E_z]

def ellipsoid_inside_check(x,y,z):
    x = x *1e8
    y = y *1e8
    z = z *1e8

    if (x**2)/(R_1**2) + (y**2)/(R_2**2) + (z**2)/(R_3**2) > 1:
        return True
    else:
        return False

def enhancement_ellipsoid(x,y,f1,L_factor_R_1,ellipsoid_permittivity):
    prefactor = (R_1*R_2*R_3/2)* ((surrounding_permittivity - ellipsoid_permittivity)/ellipsoid_permittivity)
    V = (R_1*R_2*R_3) * 1e-24
    z = 0
    r = math.sqrt(x**2 + y**2 + z**2)
    P_x = f1 / (L_factor_R_1 + 1/(ellipsoid_permittivity + surrounding_permittivity))
    if ellipsoid_inside_check(x,y,z) == False:
        E_x = -(f1 - L_factor_R_1* P_x)
        E_y = 0
    else:
        E_x = f1 * (1 + surrounding_permittivity*(ellipsoid_permittivity - surrounding_permittivity) * ((surrounding_permittivity) )/(surrounding_permittivity + L_factor_R_1 * (ellipsoid_permittivity - surrounding_permittivity)) * (V/r**5) * (2 * (x**2) - (y**2) - (z**2)))
        E_y = f1 * (surrounding_permittivity * (ellipsoid_permittivity - surrounding_permittivity) * ((surrounding_permittivity) )/(surrounding_permittivity + L_factor_R_1 * (ellipsoid_permittivity - surrounding_permittivity)) * (V/r**5) * (3*x*y))
    return math.sqrt(E_x**2/f1**2 + E_y**2/f1**2)

def oblate_definite(a,b):
    e = math.sqrt(1-b**2/a**2)
    prefactor = -1/((a**3)*e**3)
    num = 1+e 
    denom = 1-e 
    return prefactor *(2*e - math.log(num/denom))



def prolate_definite(a,c):
    x = a**2
    y = c**2
    num_1 = math.sqrt(x-y)/c
    num_2 = math.atan(math.sqrt(x-y)/c)
    denom_1 = (x-y)**(3/2)
    return 2 * (num_1 - num_2) / denom_1



if __name__ == '__main__':
    main()
