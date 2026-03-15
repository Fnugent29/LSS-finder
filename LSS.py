import sys
import numpy as np
from numpy.core.function_base import linspace
from tqdm import trange

print("Minimum mass laminate stacking sequence finder")

print("-----------------------------------------")
print("Forces")
print("-----------------------------------------")
print("")

Nx = float(input("Input the value of Nx in N/mm: "))
Ny = float(input("Input the value of Ny in N/mm: "))
Nxy = float(input("Input the value of Nxy in N/mm: "))
Mx = float(input("Input the value of Mx in Nmm/mm: "))
My = float(input("Input the value of My in Nmm/mm: "))
Mxy = float(input("Input the value of Mxy in Nmm/mm: "))
FOS = float(input("Input the Factor of Safety: "))

print("")
print("Material Properties")
print("-----------------------------------------")
print("")

E_1 = float(input("Input the value of E_1 in MPa: "))
E_2 = float(input("Input the value of E_2 in MPa: "))
G_12 = float(input("Input the value of G_12 in MPa: "))
v_12 = float(input("Input the value of v_12: "))
Xt = float(input("Input the value of Xt in MPa: "))
Xc = float(input("Input the value of Xc in MPa: "))
Yt = float(input("Input the value of Yt in MPa: "))
Yc = float(input("Input the value of Yc in MPa: "))
S = float(input("Input the value of S in MPa: "))

print("")
print("Stiffeners (X-direction)")
print("-----------------------------------------")
print("")

n_stiff_x = int(input("Input the number of stiffeners in x-direction: "))
if n_stiff_x != 0:
    n_stiff_x_type = int(input("Input stiffener type (0 for I-shape, 1 for T_shape): "))
    if n_stiff_x_type == 1:
        w1x = float(input("Input the value of w1 in mm: "))
        w2x = float(input("Input the value of w2 in mm: "))
        t1x = float(input("Input the value of t1 in mm: "))
        t2x = float(input("Input the value of t2 in mm: "))
        A1x = t1x*w1x
        A2x = t2x*w2x
    else:
        wx = float(input("Input the value of w in mm: "))
        tx = float(input("Input the value of t in mm: "))
        Ax = wx*tx

print("")
print("Stiffeners (Y-direction)")
print("-----------------------------------------")
print("")

n_stiff_y = int(input("Input the number of stiffeners in the y-direction: "))
if n_stiff_y != 0:
    n_stiff_y_type = int(input("Input the stiffener type (0 for I-shape 1 for T_shape): "))
    if n_stiff_y_type == 1:
        w1y = float(input("Input the value of w1 in mm: "))
        w2y = float(input("Input the value of w2 in mm: "))
        t1y = float(input("Input the value of t1 in mm: "))
        t2y = float(input("Input the value of t2 in mm: "))
        A1y = w1y*t1y
        A2y = w2y*t2y
    else:
        wy = float(input("Input the value of w in mm: "))
        ty = float(input("Input the value of t in mm: "))
        Ay = wy*ty

print("")
print("Stiffeners material")
print("-----------------------------------------")
print("")

if n_stiff_y != 0 or n_stiff_x != 0:
    E_stiff = float(input("Input axial stiffness of stiffener material: "))
    G_stiff = float(input("Input shear modulus of stiffener material: "))
    v_stiff = (2*G_stiff/E_stiff) - 1

print("")
print("Panel Dimensions")
print("-----------------------------------------")
print("")

panel_x = int(input("Input panel length in mm: "))
panel_y = int(input("Input panel width in mm: "))

print("")
print("Solver setup")
print("-----------------------------------------")
print("")

X = int(input("Starting ply number: "))

n_angles = 4

v_21 = E_2*v_12/E_1
Q_11 = E_1 / (1 - v_12 * v_21)
Q_12 = v_12 * E_2 / (1 - v_12 * v_21)
Q_22 = E_2 / (1 - v_12 * v_21)
Q_66 = G_12
Q_mat_local = [[Q_11, Q_12, 0],[Q_12, Q_22, 0], [0, 0, Q_66]]
ABD = np.zeros([6,6])
Forces = [Nx*FOS, Ny*FOS, Nxy*FOS, Mx*FOS, My*FOS, Mxy*FOS]
z = 0
Y = X

while Y > X - 1:
    n_plies = Y
    print("Testing laminates with ", Y, " Plies\n")
    Y += 2
    n_stacks = n_angles ** n_plies
    n_stacks_start = n_angles ** (n_plies - 2)
    thickness = 0.25
    height = n_plies * thickness
    h = linspace(-height, height, 2 * n_plies + 1)
    
    # Changing base to base number of angles
    def base_convert(a, b):
        result = np.zeros(n_plies)
        k = 0
        while a > 0:
            result[k] = a % b
            a = a // b
            k += 1
        return result


    for i in trange(n_stacks_start ,n_stacks):
        dat = base_convert(i, n_angles)
        ply_pos = 0
        ply_neg = 0
        for j in range(0, n_plies):
            if n_angles == 4:
                if dat[j] == 0:
                    ply_neg += 1
                else:
                    ply_pos += 1
                    
        # Checking that there are the same number of positive angle plies as negative angle plies
        if ply_pos == ply_neg:
            d = np.zeros(n_plies - 1)
            for j in range(1, n_plies):
                d[j - 1] = abs(dat[j - 1] - dat[j])

            k = 0
            for j in range(1, n_plies):
                if d[j - 1] > 1:
                    k += 1
                else:
                    k += 0

            L = []
            
            # Checking that consecutive plies change by no more than 45 degrees.
            if k == 0:
                d = np.zeros(n_plies - 1)
                for j in range(3, n_plies):
                    if dat[j - 3] == dat[j - 2] and dat[j - 2] == dat[j - 1] and dat[j - 1] == dat[j]:
                        d[j - 3] = 0
                    else:
                        d[j - 3] = 1

                if dat[n_plies - 3] == dat[n_plies - 2] and dat[n_plies - 2] == dat[n_plies - 1]:
                    d[n_plies - 3] = 0
                else:
                    d[n_plies - 3] = 1

                if dat[n_plies - 2] == dat[n_plies - 1]:
                    d[n_plies - 2] = 0
                else:
                    d[n_plies - 2] = 1

                k = 0
                for j in range(3, n_plies + 2):
                    if d[j - 3] == 0:
                        k += 1
                    else:
                        k += 0
                        
                # Checking that there are no more than 3 consecutive plies with the same ply angle
                if k == 0:
                    unique, counts = np.unique(dat, return_counts=True)
                    d = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0}
                    d.update(dict(zip(unique, counts)))
                    
                    # Checking all four ply angles are present
                    if len(d) > 3:
                        if n_angles == 4:
                            
                            # Checking that atleast 10% of the laminate is taken up by each ply
                            if int(d[0]) > np.floor(0.1*n_plies) and int(d[1]) > np.floor(0.1*n_plies) and int(d[2]) > np.floor(0.1*n_plies) and int(d[3]) > np.floor(0.1*n_plies):
                                if dat[0] != 1:
                                    k = 0
                                    T = np.zeros([n_plies, 3, 3])
                                    Q_mat = np.zeros([2 * n_plies, 3, 3])
                                    
                                    # Compute Transformed Q-Matrices
                                    for lam in dat:
                                        c = round(np.cos(np.pi*(lam*45 - 45)/180), 5)
                                        s = round(np.sin(np.pi*(lam*45 - 45)/180), 5)
                                        T = [[c**2, s**2, -2*c*s],[s**2, c**2, 2*c*s],[s*c, -s*c, c**2 - s**2]]
                                        T_2 = [[c**2, s**2, c*s],[s**2, c**2, -c*s],[-2*s*c, 2*s*c, c**2 - s**2]]
                                        Q_mat[k,:,:] = np.matmul(T, (np.matmul(Q_mat_local, T_2)))
                                        Q_mat[2*n_plies - k - 1, :, :] = Q_mat[k,:,:]
                                        k += 1

                                    # Compute ABD Matrix
                                    ABD = np.zeros([6, 6])
                                    for x in range(0,3):
                                        for y in range(0,3):
                                            for z in range(0, 2*n_plies):
                                                ABD[y,x] += Q_mat[z,y,x]*(h[z+1] - h[z])
                                                ABD[y,x] = round(ABD[y,x], 5)
                                                ABD[y+3,x] += (1/2)*Q_mat[z,y,x]*(h[z+1]**2 - h[z]**2)
                                                ABD[y+3,x] = round(ABD[y+3,x], 5)
                                                ABD[y,x+3] = ABD[y+3,x]
                                                ABD[y+3,x+3] += (1/3)*Q_mat[z,y,x]*(h[z+1]**3 - h[z]**3)
                                                ABD[y+3,x+3] = round(ABD[y+3,x+3], 5)
                                    
                                    # Calculate Corrections for ABD Matrix according to Smeared Stiffness Theory in X-Direction
                                    if n_stiff_x != 0:           
                                        if n_stiff_x_type == 1:
                                            
                                            #Calculating Inertias
                                            dx = ((w1x * t1x)*(w1x/2 + n_plies*thickness) + (w2x * t2x)*(w1x + t2x/2 + n_plies * thickness))/(w1x*t1x + w2x*t2x)
                                            Ixx = (w1x * t1x)*(dx - w1x/2)**2 + (t1x * w1x**3)/12 + (w2x * t2x)*(w1x + t2x/2 -dx)**2 + (t2x**3 * w2x)/12
                                            Iyy = (t1x**3 * w1x)/12 + (t2x * w2x**3)/12
                                            Ix = Ixx + (A1x + A2x)*(dx + n_plies*thickness)**2
                                            J_1 = Ixx + Iyy
                                            s = panel_x/(n_stiff_x - 1)
                                            
                                            # Ky from Timoshenko-Ehrenfest Beam Theory
                                            m = (w2x * t2x)/(w1x * t1x)
                                            n = w2x/w1x
                                            ky = (10*(1 + v_stiff)*(1+4*m)**2)/((12 + 96*m + 276*m**2 + 192*m**3) + v_stiff*(11 + 88*m + 248*m**2 + 216*m**3) + 30*n**2*(m + m**2) + 10*v_stiff*n**2*(4*m + 5*m**2 + m**3))
                                            
                                            #Updating ABD matrix
                                            ABD[0,0] += E_stiff * (A1x + A2x)/s
                                            ABD[2,2] += G_stiff * (A1x + A2x) * ky/(4*s)
                                            ABD[3,0] += E_stiff * (A1x + A2x) * (dx + n_plies*thickness)/s
                                            ABD[0,3] = ABD[3,0]
                                            ABD[2,5] += G_stiff * (A1x + A2x) * ky * (dx + n_plies*thickness)/(4*s)
                                            ABD[5,2] = ABD[2,5]
                                            ABD[3,3] += E_stiff * Ix/s
                                            ABD[5,5] += G_stiff * J_1/(4*s)
                                        else:
                                            
                                            #Calculating Inertia's
                                            Ixx = (wx**3 * tx)/12
                                            Iyy = (tx**3 * wx)/12
                                            Ix = Ixx + Ax*(wx/2)**2
                                            J_1 = Ixx + Iyy
                                            s = panel_x/(n_stiff_x - 1)
                                            
                                            # Ky from Timoshenko-Ehrenfest Beam Theory
                                            ky = 10*(1 + v_stiff)/(12 + 11*v_stiff)
                                            
                                            #Updating ABD matrix
                                            ABD[0,0] += E_stiff * Ax/s
                                            ABD[2,2] += G_stiff * Ax * ky/(4*s)
                                            ABD[3,0] += E_stiff * Ax * (wx/2 + n_plies*thickness)
                                            ABD[0,3] = ABD[3,0]
                                            ABD[2,5] += G_stiff * Ax * ky * (wx/2 + n_plies*thickness)/(4*s)
                                            ABD[5,2] = ABD[2,5]
                                            ABD[3,3] += E_stiff * Ix/s
                                            ABD[5,5] += G_stiff * J_1/s
                                            
                                    # Calculate Corrections for ABD Matrix according to Smeared Stiffness Theory in Y-Direction
                                    if n_stiff_y != 0:
                                        if n_stiff_y_type == 1:
                                            
                                            #Calculating Inertias
                                            dy = ((w1y * t1y)*(w1y/2 + n_plies*thickness) + (w2y * t2y)*(w1y + t2y/2 + n_plies * thickness))/(w1y*t1y + w2y*t2y)
                                            Ixx = (w1y * t1y)*(dy - w1y/2)**2 + (t1y * w1y**3)/12 + (w2y * t2y)*(w1y + t2y/2 -dy)**2 + (t2y**3 * w2y)/12
                                            Iyy = (t1y**3 * w1y)/12 + (t2y * w2y**3)/12
                                            Iy = Ixx + (A1y + A2y)*(dy + n_plies*thickness)**2
                                            J_2 = Ixx + Iyy
                                            s = panel_y/(n_stiff_y - 1)
                                            
                                            # Ky from Timoshenko-Ehrenfest Beam Theory
                                            m = (w2y * t2y)/(w1y * t1y)
                                            n = w2y/w1y
                                            ky = (10*(1 + v_stiff)*(1 + 4*m)**2)/((12 + 96*m + 276*m**2 + 192*m**3) + v_stiff*(11 + 88*m + 248*m**2 + 216*m**3) + 30*n**2*(m + m**2) + 10*v_stiff*n**2*(4*m + 5*m**2 + m**3))
                                            
                                            #Updating ABD matrix
                                            ABD[1,1] += E_stiff * (A1y + A2y)/s
                                            ABD[2,2] += G_stiff * (A1y + A2y) * ky/(4*s)
                                            ABD[4,1] += E_stiff * (A1y + A2y) * (dy + n_plies*thickness)/s
                                            ABD[1,4] = ABD[4,1]
                                            ABD[2,5] += G_stiff * (A1y + A2y) * ky * (dy + n_plies*thickness)/(4*s)
                                            ABD[5,2] = ABD[2,5]
                                            ABD[4,4] += E_stiff * Iy/s
                                            ABD[5,5] += G_stiff * J_2/s
                                        else:
                                            
                                            #Calculating Inertias
                                            Ixx = (wy**3 * ty)/12
                                            Iyy = (ty**3 * wy)/12
                                            Iy = Ixx + Ay*(wy/2 + n_plies*thickness)
                                            J_2 = Ixx + Iyy
                                            s = panel_y/(n_stiff_y - 1)
                                            
                                            # Ky from Timoshenko-Ehrenfest Beam Theory
                                            ky = 10*(1 + v_stiff)/(12 + 11*v_stiff)
                                            
                                            #Updating ABD matrix
                                            ABD[1,1] += E_stiff * Ay/s
                                            ABD[2,2] += G_stiff * Ay * ky/(4*s)
                                            ABD[4,1] += E_stiff * Ay * (wy/2 + n_plies*thickness)
                                            ABD[1,4] = ABD[4,1]
                                            ABD[2,5] += G_stiff * Ay * ky * (wx/2 + n_plies*thickness)/(4*s)
                                            ABD[5,2] = ABD[2,5]
                                            ABD[4,4] += E_stiff * Iy/s
                                            ABD[5,5] += G_stiff * J_2/s
                                            
                                    #Compute Strains and Curvatures
                                    strain_curvature = np.matmul(Forces, np.linalg.inv(ABD))
                                    k = 0
                                    hoffman_1 = 1
                                    hoffman_2 = 1
                                    hoffman_3 = 1
                                    hoffman_4 = 1
                                    strain_mat_x = np.zeros(4*n_plies)
                                    strain_mat_y = np.zeros(4*n_plies)
                                    strain_mat_xy = np.zeros(4*n_plies)
                                    principle_stress_mat_x = np.zeros(4*n_plies)
                                    principle_stress_mat_y = np.zeros(4*n_plies)
                                    principle_stress_mat_xy = np.zeros(4*n_plies)
                                    hoffman_mat = np.zeros(4*n_plies)
                                    
                                    for lam in dat:
                                        if k == 0 or (hoffman_1 <= 1 and hoffman_2 <= 1 and hoffman_3 <= 1 and hoffman_4 <= 1):
                                            
                                            # Local Strain Calculation
                                            strain_1 = strain_curvature[0:3] + strain_curvature[3:6] * h[k]
                                            strain_2 = strain_curvature[0:3] + strain_curvature[3:6] * h[k + 1]
                                            strain_3 = strain_curvature[0:3] + strain_curvature[3:6] * -h[k]
                                            strain_4 = strain_curvature[0:3] + strain_curvature[3:6] * -h[k + 1]
                                            
                                            #Rotation Matrices
                                            c = round(np.cos(np.pi*(lam*45 - 45) / 180), 5)
                                            s = round(np.sin(np.pi*(lam*45 - 45) / 180), 5)
                                            T_3 = [[c ** 2, s ** 2, 2 * c * s], [s ** 2, c ** 2, -2 * c * s],[-s * c, s * c, c ** 2 - s ** 2]]
                                            T = [[c ** 2, s ** 2, -2 * c * s], [s ** 2, c ** 2, 2 * c * s],
                                                 [s * c, -s * c, c ** 2 - s ** 2]]
                                            T_2 = [[c ** 2, s ** 2, c * s], [s ** 2, c ** 2, -c * s],
                                                   [-2 * s * c, 2 * s * c, c ** 2 - s ** 2]]
                                                   
                                            # Transformed Q-Matrix
                                            Q_mat = np.matmul(T, (np.matmul(Q_mat_local, T_2)))
                                            
                                            # Global Stresses
                                            stress_1 = np.matmul(Q_mat, strain_1)
                                            stress_2 = np.matmul(Q_mat, strain_2)
                                            stress_3 = np.matmul(Q_mat, strain_3)
                                            stress_4 = np.matmul(Q_mat, strain_4)
                                            
                                            #Principal Stresses
                                            principle_stress_mat_x[2*k] = np.matmul(T_3, stress_1)[0]
                                            principle_stress_mat_y[2*k] = np.matmul(T_3, stress_1)[1]
                                            principle_stress_mat_xy[2*k] = np.matmul(T_3, stress_1)[2]
                                            principle_stress_1 = np.matmul(T_3, stress_1)
                                            principle_stress_mat_x[2*k+1] = np.matmul(T_3, stress_2)[0]
                                            principle_stress_mat_y[2*k+1] = np.matmul(T_3, stress_2)[1]
                                            principle_stress_mat_xy[2*k+1] = np.matmul(T_3, stress_2)[2]
                                            principle_stress_2 = np.matmul(T_3, stress_2)
                                            principle_stress_mat_x[4*n_plies - 2*k - 1] = np.matmul(T_3, stress_3)[0]
                                            principle_stress_mat_y[4*n_plies - 2*k - 1] = np.matmul(T_3, stress_3)[1]
                                            principle_stress_mat_xy[4*n_plies - 2*k - 1] = np.matmul(T_3, stress_3)[2]
                                            principle_stress_3 = np.matmul(T_3, stress_3)
                                            principle_stress_mat_x[4*n_plies - 2*k - 2] = np.matmul(T_3, stress_4)[0]
                                            principle_stress_mat_y[4*n_plies - 2*k - 2] = np.matmul(T_3, stress_4)[1]
                                            principle_stress_mat_xy[4*n_plies - 2*k - 2] = np.matmul(T_3, stress_4)[2]
                                            principle_stress_4 = np.matmul(T_3, stress_4)
                                            
                                            # Hoffman Failure Criterion
                                            hoffman_1 = (1 / Xt - 1 / Xc) * principle_stress_1[0] + (1 / Yt - 1 / Yc) * \
                                                      principle_stress_1[1] + principle_stress_1[0] ** 2 / (Xt * Xc) + principle_stress_1[
                                                          1] ** 2 / (Yt * Yc) - (principle_stress_1[0] * principle_stress_2[1]) / (
                                                                  Xt * Xc) + (principle_stress_1[2] / S) ** 2
                                            hoffman_2 = (1 / Xt - 1 / Xc) * principle_stress_2[0] + (1 / Yt - 1 / Yc) * \
                                                      principle_stress_2[1] + principle_stress_2[0] ** 2 / (Xt * Xc) + principle_stress_2[
                                                          1] ** 2 / (Yt * Yc) - (principle_stress_2[0] * principle_stress_2[1]) / (
                                                                  Xt * Xc) + (principle_stress_2[2] / S) ** 2
                                            hoffman_3 = (1 / Xt - 1 / Xc) * principle_stress_3[0] + (1 / Yt - 1 / Yc) * \
                                                      principle_stress_3[1] + principle_stress_3[0] ** 2 / (Xt * Xc) + principle_stress_3[
                                                          1] ** 2 / (Yt * Yc) - (principle_stress_3[0] * principle_stress_2[1]) / (
                                                                  Xt * Xc) + (principle_stress_3[2] / S) ** 2
                                            hoffman_4 = (1 / Xt - 1 / Xc) * principle_stress_4[0] + (1 / Yt - 1 / Yc) * \
                                                      principle_stress_4[1] + principle_stress_4[0] ** 2 / (Xt * Xc) + principle_stress_4[
                                                          1] ** 2 / (Yt * Yc) - (principle_stress_4[0] * principle_stress_4[1]) / (
                                                                  Xt * Xc) + (principle_stress_4[2] / S) ** 2
                                            hoffman_mat[2*k] = hoffman_1
                                            hoffman_mat[2*k +1] = hoffman_2
                                            hoffman_mat[4*n_plies - 2*k - 1] = hoffman_3
                                            hoffman_mat[4*n_plies - 2*k - 2] = hoffman_4
                                            k += 1
                                            
                                    # If all plies don't exceed' the failure criterion'
                                    if k == n_plies:
                                        D_11 = ABD[3,3]
                                        D_12 = ABD[3,4]
                                        D_22 = ABD[4,4]
                                        D_66 = ABD[5,5]
                                        
                                        N_x_cr = 0
                                        N_y_cr = 0
                                        
                                        # Performing Buckling Analysis
                                        for m in range(1, 20):
                                            for n in range(1, 20):
                                                M = (np.pi/(m**2 * panel_x**2))*(D_11*m**4 + 2*(D_12 + 2*D_66)*m**2 * n**2 *(panel_x/panel_y)**2 + D_22*n**4 * (panel_x/panel_y)**4)
                                                if N_x_cr > M or (n == 1 and m == 1):
                                                    N_x_cr = M
                                                
                                        # Calculating critical buckling forces in other directions
                                        N_y_cr = N_x_cr*Forces[1]/Forces[0]
                                        N_xy_cr = N_x_cr * Forces[2] / Forces[0]
                                        
                                        # Checking critical buckling forces are below the applied loads
                                        if Forces[0] < N_x_cr and Forces[1] < N_y_cr and Forces[2] < N_xy_cr:
                                        
                                            print("Lay-up: ", dat, ABD[3:6, 3:6], " Has PASSED")
                                            print("---------------------------------------------------------------")
                                            print("")
                                            print("Hoffman Criterion: ", hoffman_mat)
                                            print("---------------------------------------------------------------")
                                            print("")
                                            print("Principal Stress 1: ", principle_stress_mat_x)
                                            print("Principal Stress 2: ", principle_stress_mat_y)
                                            print("Shear Stress: ", principle_stress_mat_xy)
                                            print("---------------------------------------------------------------")
                                            print("")
                                            print("Critical Buckling X: ", N_x_cr)
                                            print("Critical Buckling Y: ", N_y_cr)
                                            print("Critical Buckling XY: ", N_xy_cr)
                                        
                                            response = int(input("Is this a good laminate? (1 for Yes, 0 for No): "))
                                            if response == 1:
                                                sys.exit("Lay-up has been found")


