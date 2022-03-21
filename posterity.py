def prolate_indefinite(xi,a,b):
    x = a**2
    y = b**2
    mp.dps = 25
    mp.pretty = True
    #note only defnined for abs (xi + y)/(y-x) < 1
    num_1 = 2 * hyp2f1(-0.5,1,0.5,(xi + y)/(y-x))
    denom_1 = math.sqrt(xi + y) * (x - y)
    #print(num_1/denom_1)
    return num_1/denom_1


def oblate_indefinite(xi,a,b):
    x = a**2
    y = b**2
    num_1 = math.sqrt(xi + y)
    num_2 = math.atan(math.sqrt(xi + y)/math.sqrt(x-y))
    denom_1 = (xi + x) * (xi - y)
    denom_2 = (x-y)**(3/2)
    #print((num_1/denom_1) + (num_2/denom_2))
    return (num_1/denom_1) + (num_2/denom_2)

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

#Returns the electric field for the ellispoidal particle 
def E_ellipsoid_inside(L_factor_R_1, L_factor_R_2,L_factor_R_3,incident_field):
    ellipsoid_field = np.zeros((timesteps, 3))
    for t in range(timesteps -1):
        current_field = incident_field[t]
        P_x = e_0 * (sphere_permittivity - surrounding_permittivity) * ((surrounding_permittivity) * current_field[0])/(surrounding_permittivity + L_factor_R_1 (sphere_permittivity - surrounding_permittivity))
        P_y = e_0 * (sphere_permittivity - surrounding_permittivity) * ((surrounding_permittivity) * current_field[1])/(surrounding_permittivity + L_factor_R_1 (sphere_permittivity - surrounding_permittivity))
        P_z = e_0 * (sphere_permittivity - surrounding_permittivity) * ((surrounding_permittivity) * current_field[2])/(surrounding_permittivity + L_factor_R_1 (sphere_permittivity - surrounding_permittivity))
        ellipsoid_field[t][0] = current_field[0] - L_factor_R_1* P_x/(e_0 * surrounding_permittivity)
        ellipsoid_field[t][1] = current_field[1] - L_factor_R_2* P_y/(e_0 * surrounding_permittivity)
        ellipsoid_field[t][2] = current_field[2] - L_factor_R_3* P_z/(e_0 * surrounding_permittivity)
    return ellipsoid_field

def Xi(x,y,z):
    x = x *1e9
    y = y *1e9
    z = z *1e9
    h = math.sqrt(R_1**2 - R_2**2)
    k = math.sqrt(R_1**2 - R_3**2)
    a_1 = -((x**2)+(y**2)+(z**2)+(h**2)+(k**2))
    a_2 = (x**2)*((h**2) + (k**2)) + (y**2)*(k**2) + (z**2)*(h**2) + (h**2)* (k**2)
    a_3 = (x**2)*(h**2)*(k**2)
    Q = ((a_1**2) - 3 * a_2) / 9
    R = ((9*a_1*a_2) - (27 * a_3) - (2 * a_1**3))/54
    corrected_theta = R/math.sqrt(Q**3)
    if abs(corrected_theta)>1:
        corrected_theta = corrected_theta/ abs(corrected_theta)
    #print(corrected_theta)
    theta = math.acos(corrected_theta)
    xi = 2 * math.sqrt(Q) * math.cos (theta/3) - (a_1/3)
    return math.sqrt(xi)