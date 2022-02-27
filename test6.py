import math
import numpy as np
import random as rand
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import os
from matplotlib import animation
#How the enhancement varies with different resonance conditions
#Allow for selection from an angle from the x axis
#plot the x and y components of the colourplot
#add colourbars
#run for more cycles

#Input parameters in SI units
#sphere radius 50nm
sphere_radius = 50e-9
#wavelength 500nm
wavelength = 500e-9
#Amplitude V/m
#1mw comfortable level -> 1 micron assuming some power density -> known flux
amplitude = 1e5
sphere_permittivity = -2.0001
surrounding_permittivity = 1
c = 299792458 
e = 1.60217662e-19
direction = [1,1,0]
timesteps = 200
electron_mass = 9.10938356e-31
omega = 2 * math.pi * c/ wavelength
#setting phi and theta and r
# theta [0,pi] phi [0,2pi]
phi = 0
theta = 0
#0.1 nm away from the sphere radius
r = sphere_radius + 0.0000000001
phi = phi%math.pi
theta = theta%(2*math.pi)



def main():
    timestep = interpolate()
    time = 0
    incident_field = np.zeros((timesteps,3))
    #stores the incident field in an array
    for t in range(timesteps):
        Ex = amplitude* math.cos(omega*time) * direction[0]
        Ey = amplitude* math.cos(omega*time) * direction[1]
        Ez = amplitude* math.cos(omega*time) * direction[2]
        incident_field[t][0] = Ex
        incident_field[t][1] = Ex
        incident_field[t][2] = Ex
        time += timestep
    #generating coordinate points for electron
    electron_position = np.zeros((timesteps,3))
    electron_velocity = np.zeros((timesteps,3))
    x = r * math.cos(theta)
    y = r * math.cos(phi) * math.sin(theta)
    z = r * math.sin(phi) * math.sin(theta)
    electron_position[0][0] = x
    electron_position[0][1] = y
    electron_position[0][2] = z
    graph_field = []
    #calculating the electric field at the electron position and updating the position of the electron
    for t in range(timesteps -1):
        current_field = incident_field[t]
        field = electric_field(electron_position[t][0],electron_position[t][1],electron_position[t][2],current_field)
        graph_field.append(field[0])
        if t == 0: 
            electron_position[1][0] = electron_position[0][0] + (field[0]* e/electron_mass) * (timestep**2)
            electron_position[1][1] = electron_position[0][1] + (field[1]* e/electron_mass) * (timestep**2)
            electron_position[1][2] = electron_position[0][2] + (field[2]* e/electron_mass) * (timestep**2)
        else:
            electron_position[t+1][0] = (2 * electron_position[t][0]) - electron_position[t-1][0] + (timestep**2)  * (field[0] * e/electron_mass)
            electron_position[t+1][1] = (2 * electron_position[t][1]) - electron_position[t-1][1] + (timestep**2)  * (field[1] * e/electron_mass)
            electron_position[t+1][2] = (2 * electron_position[t][2]) - electron_position[t-1][2] + (timestep**2)  * (field[2] * e/electron_mass)
    #calculating the velocity using the verlet algorithm
    for t in range(timesteps -1):
        if t == 0: 
            electron_velocity[0][0] = 0
            electron_velocity[0][1] = 0
            electron_velocity[0][2] = 0
        else: 
            electron_velocity[t][0] = (electron_position[t+1][0] - electron_position[t-1][0]) / (2 * timestep)
            electron_velocity[t][1] = (electron_position[t+1][1] - electron_position[t-1][1]) / (2 * timestep)
            electron_velocity[t][2] = (electron_position[t+1][2] - electron_position[t-1][2]) / (2 * timestep)

    timestep_graph(electron_position, timestep,graph_field, electron_velocity)
    #Creating a grid of x,y plots
    E2= np.vectorize(E)
    nx,ny= 100,100 #level of discretisation
    xr = np.linspace(-5*sphere_radius,5*sphere_radius,nx)
    yr = np.linspace(-5*sphere_radius,5*sphere_radius,ny)
    X,Y = np.meshgrid(xr,yr)
    plot_Ex = np.zeros((nx,ny))
    plot_Ey = np.zeros((nx,ny))
    f1 = incident_field[t][0]
    ex,ey = E2(x=X, y = Y,f1 = f1) 
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
    Q = ax.quiver(xr, yr, plot_Ex, plot_Ey, pivot = 'mid')
    #plt.show()
    #add the nanoparticle to the image
    #ax.add_artist(Circle((0,0),sphere_radius))
    anim = animation.FuncAnimation(fig, update_quiver, fargs=(Q, X, Y, incident_field, E2),
                               interval=500, blit=False)
    fig.tight_layout()
    # saving to mp4 using ffmpeg writer
    anim.save('Electric field.gif', fps = 60, dpi = 300)
    plt.close()
    


def inside_check(x,y,z):
    r = math.sqrt((x**2) + (y**2) + (z**2))
    if r < sphere_radius: 
        return True
    else:
        return False
#Returns the electric field for the quiver plot
def E(x,y,f1):
    z = 0
    r = math.sqrt((x**2) + (y**2))
    
    if inside_check(x,y,z) == False: 
        Ex = (1 + ((sphere_permittivity - surrounding_permittivity)/(sphere_permittivity + 2* surrounding_permittivity)) * ((sphere_radius**3)/r**5) * (2* (x**2) -(y**2) - (z**2)))* f1
        Ey = ((sphere_permittivity - surrounding_permittivity)/(sphere_permittivity + 2* surrounding_permittivity)) * ((sphere_radius**3)/r**5) * (3*x*y) * f1
        return Ex,Ey
    if inside_check(x,y,z) == True:
        Ex = f1 * (3*surrounding_permittivity/(sphere_permittivity + 2* surrounding_permittivity))
        Ey = 0
        return Ex,Ey

#Returns the electric field for the electron
def electric_field(x,y,z,current_field):
    r = math.sqrt((x**2) + (y**2) + (z**2))
    if inside_check(x,y,z) == False: 
        Ex = (1 + ((sphere_permittivity - surrounding_permittivity)/(sphere_permittivity + 2* surrounding_permittivity))
         * ((sphere_radius**3)/r**5) * (2* (x**2) -(y**2) - (z**2)))* current_field[0]
        Ey = ((sphere_permittivity - surrounding_permittivity)/(sphere_permittivity + 2* surrounding_permittivity)) * ((sphere_radius**3)/r**5) * (3*x*y) * current_field[0]
        Ez = ((sphere_permittivity - surrounding_permittivity)/(sphere_permittivity + 2* surrounding_permittivity)) * ((sphere_radius**3)/r**5) * (3*x*z) * current_field[0]
        return[Ex,Ey,Ez]
    if inside_check(x,y,z) == True:
        Ex = current_field[0] * (3*surrounding_permittivity/(sphere_permittivity + 2* surrounding_permittivity))
        Ey = 0
        Ez = 0
        return[Ex,Ey,Ez]
#function that plots different graphs
def timestep_graph(electron_position,timestep,graph_field,electron_velocity):
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
    plt.plot(timestep_for_graph,xpos)
    plt.savefig("Xposition of electron.png", dpi = 300, bbox_inches = 'tight')
    plt.close()
    plt.figure(1)
    plt.plot(timestep_for_graph,ypos)
    plt.savefig("Yposition of electron.png", dpi = 300, bbox_inches = 'tight')
    plt.close()
    plt.figure(2)
    plt.plot(timestep_for_graph,zpos)
    plt.savefig("Zposition of electron.png", dpi = 300, bbox_inches = 'tight')
    plt.close()
    plt.figure(3)
    plt.plot(timestep_for_graph,graph_field)
    plt.savefig("electric field amplitude.png", dpi = 300, bbox_inches = 'tight')
    plt.close()
    plt.plot(timestep_for_graph,xvel)
    plt.savefig("Xvelocity of electron.png", dpi = 300, bbox_inches = 'tight')
    plt.close()
    plt.plot(timestep_for_graph,yvel)
    plt.savefig("Yvelocity of electron.png", dpi = 300, bbox_inches = 'tight')
    plt.close()
    plt.plot(timestep_for_graph,zvel)
    plt.savefig("Zvelocity of electron.png", dpi = 300, bbox_inches = 'tight')
    plt.close()

#function calculates the timestep
def interpolate():
    hertz = c/wavelength
    timestep = 1/(hertz * 10)
    return timestep

#quiverplot animation
def update_quiver(num,Q,X,Y,incident_field,E2):
    f1 = incident_field[num][0]
    ex,ey = E2(x=X, y = Y,f1 = f1)
    Q.set_UVC(ex,ey)
    return Q


if __name__ == '__main__':
    main()