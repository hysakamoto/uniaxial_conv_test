import numpy as np
from shape_functions import *
from math import sqrt
import sys

## theyt represent the orientations of each branch
xmax_sign = [-1,  0, 0,-1,-1, 0, 0,-1]
xmin_sign = [ 0, +1,+1, 0, 0,+1,+1, 0]
ymax_sign = [-1, -1, 0, 0,-1,-1, 0, 0]
ymin_sign = [ 0,  0,+1,+1, 0, 0,+1,+1]
zmax_sign = [-1, -1,-1,-1, 0, 0, 0, 0]
zmin_sign = [ 0,  0, 0, 0,+1,+1,+1,+1]

## hex8 numerical integration points
p_int = [[+sqrt(1/3), -sqrt(1/3), -sqrt(1/3)], \
         [+sqrt(1/3), +sqrt(1/3), -sqrt(1/3)], \
         [-sqrt(1/3), +sqrt(1/3), -sqrt(1/3)], \
         [-sqrt(1/3), -sqrt(1/3), -sqrt(1/3)], \
         [+sqrt(1/3), -sqrt(1/3), +sqrt(1/3)], \
         [+sqrt(1/3), +sqrt(1/3), +sqrt(1/3)], \
         [-sqrt(1/3), +sqrt(1/3), +sqrt(1/3)], \
         [-sqrt(1/3), -sqrt(1/3), +sqrt(1/3)] ];
weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]


def integrate_hex8( N0, N1, N2, N3, N4, N5, N6, N7, E ):
    "integrate some vlaue E  over the hex8 element"
    # N0~N7: element nodes
    # E: nodal values of E

    for i in range (8):
        E[i] = E[i]*E[i]

    E_int = [ interpolate(p_int[0], E), \
              interpolate(p_int[1], E), \
              interpolate(p_int[2], E), \
              interpolate(p_int[3], E), \
              interpolate(p_int[4], E), \
              interpolate(p_int[5], E), \
              interpolate(p_int[6], E), \
              interpolate(p_int[7], E)]

    Ns = np.array([N0, N1, N2, N3, N4, N5, N6, N7])

    Value = 0
    for i in range(8):
        J = Jacobian (p_int[i][0],p_int[i][1],p_int[i][2], Ns)
        detJ = np.linalg.det(J)

        Value = Value + \
                E_int[i] * weights[i] * detJ;

    return Value
    


def integrate_hex8_xyzp( Ns, Es ):
    "integrate some vlaues Es over the hex8 element"
    # Ns: element nodes
    # E: nodal values of E

    for nv in range(4):
        for i in range (8):
            Es[nv][i] = Es[nv][i]*Es[nv][i]

    E_int = [0.0]*4
    for nv in range(4):
        E_int[nv] = [ interpolate(p_int[0], Es[nv]), \
                      interpolate(p_int[1], Es[nv]), \
                      interpolate(p_int[2], Es[nv]), \
                      interpolate(p_int[3], Es[nv]), \
                      interpolate(p_int[4], Es[nv]), \
                      interpolate(p_int[5], Es[nv]), \
                      interpolate(p_int[6], Es[nv]), \
                      interpolate(p_int[7], Es[nv])]

    Value = 4*[0.0]
    for i in range(8):
        J = Jacobian (p_int[i][0],p_int[i][1],p_int[i][2], np.array(Ns))
        detJ = np.linalg.det(J)

        for nv in range(4):
            Value[nv] = Value[nv] + \
                    E_int[nv][i] * weights[i] * detJ;

    return Value
    



## Space integration of the error over the same mesh.
def space_integrate_error( xyzp_1, xyzp_2, nodes, elems ):
    "integrate the error between xyzp_1 (reference) and xyzp_2."

    Ex_error = 0.0
    Ey_error = 0.0
    Ez_error = 0.0
    Ep_error = 0.0
    
    pr_c = 0.0
    for i in range(len(elems)):
        # show progress
        # print i, len(elems)/10.0*pr_c
        if i > len(elems)/10.0*pr_c:
            pr_c = pr_c+1
            sys.stdout.write('*')

        Ex = [xyzp_1[elems[i][0]][0] - xyzp_2[elems[i][0]][0], \
              xyzp_1[elems[i][1]][0] - xyzp_2[elems[i][1]][0], \
              xyzp_1[elems[i][2]][0] - xyzp_2[elems[i][2]][0], \
              xyzp_1[elems[i][3]][0] - xyzp_2[elems[i][3]][0], \
              xyzp_1[elems[i][4]][0] - xyzp_2[elems[i][4]][0], \
              xyzp_1[elems[i][5]][0] - xyzp_2[elems[i][5]][0], \
              xyzp_1[elems[i][6]][0] - xyzp_2[elems[i][6]][0], \
              xyzp_1[elems[i][7]][0] - xyzp_2[elems[i][7]][0]]

        Ey = [xyzp_1[elems[i][0]][1] - xyzp_2[elems[i][0]][1], \
              xyzp_1[elems[i][1]][1] - xyzp_2[elems[i][1]][1], \
              xyzp_1[elems[i][2]][1] - xyzp_2[elems[i][2]][1], \
              xyzp_1[elems[i][3]][1] - xyzp_2[elems[i][3]][1], \
              xyzp_1[elems[i][4]][1] - xyzp_2[elems[i][4]][1], \
              xyzp_1[elems[i][5]][1] - xyzp_2[elems[i][5]][1], \
              xyzp_1[elems[i][6]][1] - xyzp_2[elems[i][6]][1], \
              xyzp_1[elems[i][7]][1] - xyzp_2[elems[i][7]][1]]

        Ez = [xyzp_1[elems[i][0]][2] - xyzp_2[elems[i][0]][2], \
              xyzp_1[elems[i][1]][2] - xyzp_2[elems[i][1]][2], \
              xyzp_1[elems[i][2]][2] - xyzp_2[elems[i][2]][2], \
              xyzp_1[elems[i][3]][2] - xyzp_2[elems[i][3]][2], \
              xyzp_1[elems[i][4]][2] - xyzp_2[elems[i][4]][2], \
              xyzp_1[elems[i][5]][2] - xyzp_2[elems[i][5]][2], \
              xyzp_1[elems[i][6]][2] - xyzp_2[elems[i][6]][2], \
              xyzp_1[elems[i][7]][2] - xyzp_2[elems[i][7]][2]]

        Ep = [xyzp_1[elems[i][0]][3] - xyzp_2[elems[i][0]][3], \
              xyzp_1[elems[i][1]][3] - xyzp_2[elems[i][1]][3], \
              xyzp_1[elems[i][2]][3] - xyzp_2[elems[i][2]][3], \
              xyzp_1[elems[i][3]][3] - xyzp_2[elems[i][3]][3], \
              xyzp_1[elems[i][4]][3] - xyzp_2[elems[i][4]][3], \
              xyzp_1[elems[i][5]][3] - xyzp_2[elems[i][5]][3], \
              xyzp_1[elems[i][6]][3] - xyzp_2[elems[i][6]][3], \
              xyzp_1[elems[i][7]][3] - xyzp_2[elems[i][7]][3]]

        Es_erros = integrate_hex8_xyzp( [nodes[elems[i][0]], \
                                         nodes[elems[i][1]], \
                                         nodes[elems[i][2]], \
                                         nodes[elems[i][3]], \
                                         nodes[elems[i][4]], \
                                         nodes[elems[i][5]], \
                                         nodes[elems[i][6]], \
                                         nodes[elems[i][7]]], \
                                        [Ex, Ey, Ez, Ep] )

        Ex_error = Ex_error + Es_erros[0]
        Ey_error = Ey_error + Es_erros[1]
        Ez_error = Ez_error + Es_erros[2]
        Ep_error = Ep_error + Es_erros[3]

    Ex_error = sqrt(Ex_error)
    Ey_error = sqrt(Ey_error)
    Ez_error = sqrt(Ez_error)
    Ep_error = sqrt(Ep_error)

    return (Ex_error, Ey_error, Ez_error, Ep_error)


## Calculate the point error.
def point_error( xyzp_1, xyzp_2, nodes, elems ):
    "calculate the point error between xyzp_1 (reference) and xyzp_2."

    Ex_error = 0.0
    Ey_error = 0.0
    Ez_error = 0.0
    Ep_error = 0.0

    # 4: the initial node at (0,0,1)
    Ex_error = abs(xyzp_1[6][0] - xyzp_2[6][0])
    Ey_error = abs(xyzp_1[6][1] - xyzp_2[6][1])
    # 4: the initial node at (0,0,1)
    Ez_error = abs(xyzp_1[5][2] - xyzp_2[5][2])
    Ep_error = abs(xyzp_1[4][3] - xyzp_2[4][3])
            
    return (Ex_error, Ey_error, Ez_error, Ep_error)


## build an octree for the point -> element map
def build_octree( elem_num, \
                  branch, \
                  elems, \
                  nodes, \
                  xyz_range, \
                  xlen, ylen, zlen, \
                  depth ):

    if branch == []:
       branch = [[], [], [], [], [], [], [], []]

    elem_dist = [[], [], [], [], [], [], [], []]

    hxlen = xlen/2.0
    hylen = ylen/2.0
    hzlen = zlen/2.0

    # print "depth", depth
    # print "range", [xmin, xmax, ymin, ymax, zmin, zmax]

    for nn in range(8):

        nx = nodes[elems[elem_num][nn]][0]
        ny = nodes[elems[elem_num][nn]][1]
        nz = nodes[elems[elem_num][nn]][2]

        # assign the element to corresponding branch
        if (nx<=xyz_range[1]) & (nx >= (xyz_range[0]+hxlen)):
            if (ny >= xyz_range[2]) & (ny<=xyz_range[2]+hylen):
                if (nz >= xyz_range[4]) & (nz<=xyz_range[4]+hzlen):
                    elem_dist[1].append(elem_num)
                if (nz <= xyz_range[5]) & (nz >= (xyz_range[4]+hzlen)):
                    elem_dist[5].append(elem_num)
            if  (ny<=xyz_range[3]) & (ny >= (xyz_range[2]+hylen)):
                if (nz >= xyz_range[4]) & (nz <= (xyz_range[4]+hzlen)):
                    elem_dist[2].append(elem_num)
                if (nz<=xyz_range[5]) & (nz >= (xyz_range[4]+hzlen)):
                    elem_dist[6].append(elem_num)
        if (nx >= xyz_range[0]) & (nx <= xyz_range[0]+hxlen):
            if (ny >= xyz_range[2]) & (ny<=xyz_range[2]+hylen):
                if (nz >= xyz_range[4]) & (nz<=(xyz_range[4]+hzlen)):
                    elem_dist[0].append(elem_num) 
                if (nz >= (xyz_range[4]+hzlen)) & (nz<=xyz_range[5]):
                    elem_dist[4].append(elem_num)
            if (ny<=xyz_range[3]) & (ny >= (xyz_range[2]+hylen)):
                if (nz >= xyz_range[4]) & (nz<=xyz_range[4]+hzlen):
                    elem_dist[3].append(elem_num)
                if (nz<=xyz_range[5]) & (nz >= (xyz_range[4]+hzlen)):
                    elem_dist[7].append(elem_num)

    ## termination criteria
    if (depth-1) <= 0:
        for nn in range(8):
            if elem_dist[nn] != []:
                branch[nn].append(elem_num)
        return branch

    else:
        for dn in range(8):
            if elem_dist[dn] != []:
                branch[dn] = \
                build_octree( elem_num, \
                              branch[dn], \
                              elems, \
                              nodes, \
                              [xyz_range[0]+(xmin_sign[dn])*hxlen, \
                               xyz_range[1]+(xmax_sign[dn])*hxlen, \
                               xyz_range[2]+(ymin_sign[dn])*hylen, \
                               xyz_range[3]+(ymax_sign[dn])*hylen, \
                               xyz_range[4]+(zmin_sign[dn])*hzlen, \
                               xyz_range[5]+(zmax_sign[dn])*hzlen ], \
                              hxlen, hylen, hzlen, \
                              depth-1 )
            
        return branch


def search_branch( px, py, pz,  \
                   xyz_range, \
                   xlen, ylen, zlen, depth):

    hxlen = xlen/2.0
    hylen = ylen/2.0
    hzlen = zlen/2.0

    branch_num = -1
    
    # search the node in each branch
    if (px<=xyz_range[1]) & (px >= (xyz_range[0]+hxlen)):
        if (py >= xyz_range[2]) & (py<=xyz_range[2]+hylen):
            if (pz >= xyz_range[4]) & (pz<=xyz_range[4]+hzlen):
                branch_num = 1
            elif (pz <= xyz_range[5]) & (pz >= (xyz_range[4]+hzlen)):
                branch_num = 5
        elif  (py<=xyz_range[3]) & (py >= (xyz_range[2]+hylen)):
            if (pz >= xyz_range[4]) & (pz <= (xyz_range[4]+hzlen)):
                branch_num = 2
            elif (pz<=xyz_range[5]) & (pz >= (xyz_range[4]+hzlen)):
                branch_num = 6
    elif (px >= xyz_range[0]) & (px <= xyz_range[0]+hxlen):
        if (py >= xyz_range[2]) & (py<=xyz_range[2]+hylen):
            if (pz >= xyz_range[4]) & (pz<=(xyz_range[4]+hzlen)):
                branch_num = 0
            elif (pz >= (xyz_range[4]+hzlen)) & (pz<=xyz_range[5]):
                branch_num = 4
        elif (py<=xyz_range[3]) & (py >= (xyz_range[2]+hylen)):
            if (pz >= xyz_range[4]) & (pz<=xyz_range[4]+hzlen):
                branch_num = 3
            elif (pz<=xyz_range[5]) & (pz >= (xyz_range[4]+hzlen)):
                branch_num = 7

    if branch_num == -1:
        print "error, the point cannot be found in the range"
    elif depth <= 1:
        return [branch_num]
    else:
        branch_chain = search_branch( px, py, pz, \
            [xyz_range[0]+(xmin_sign[branch_num])*hxlen, \
             xyz_range[1]+(xmax_sign[branch_num])*hxlen, \
             xyz_range[2]+(ymin_sign[branch_num])*hylen, \
             xyz_range[3]+(ymax_sign[branch_num])*hylen, \
             xyz_range[4]+(zmin_sign[branch_num])*hzlen, \
             xyz_range[5]+(zmax_sign[branch_num])*hzlen ] , \
            hxlen, hylen, hzlen, \
            depth-1)
        branch_chain.append(branch_num)
        return branch_chain



def interpolate( xis, E ):
    "interpolate the value E at master coordinates xis"
    # E: nodeal values
    # xis: master coordinates

    PHI = [PHI0(xis[0],xis[1],xis[2]), \
           PHI1(xis[0],xis[1],xis[2]), \
           PHI2(xis[0],xis[1],xis[2]), \
           PHI3(xis[0],xis[1],xis[2]), \
           PHI4(xis[0],xis[1],xis[2]), \
           PHI5(xis[0],xis[1],xis[2]), \
           PHI6(xis[0],xis[1],xis[2]), \
           PHI7(xis[0],xis[1],xis[2])]
    
    Value = 0.0
    for i in range(8):
        Value = Value + E[i] * PHI[i]

    return Value


def Residual( xi,eta,mu, Ns, P):

    PHI = np.array([[PHI0(xi,eta,mu)], \
                 [PHI1(xi,eta,mu)], \
                 [PHI2(xi,eta,mu)], \
                 [PHI3(xi,eta,mu)], \
                 [PHI4(xi,eta,mu)], \
                 [PHI5(xi,eta,mu)], \
                 [PHI6(xi,eta,mu)], \
                 [PHI7(xi,eta,mu)]])

    # print PHI
    # print dot(Ns[0], PHI)
    # print dot(Ns[1], PHI)
    # print dot(Ns[2], PHI)

    Ns = Ns.T
    return np.array([P[0] - np.dot(Ns[0], PHI), \
                     P[1] - np.dot(Ns[1], PHI), \
                     P[2] - np.dot(Ns[2], PHI)]) 


def Jacobian (xi,eta,mu, Ns):

    dPHIs = np.array([dPHI0(xi,eta,mu),\
                      dPHI1(xi,eta,mu),\
                      dPHI2(xi,eta,mu),\
                      dPHI3(xi,eta,mu),\
                      dPHI4(xi,eta,mu),\
                      dPHI5(xi,eta,mu),\
                      dPHI6(xi,eta,mu),\
                      dPHI7(xi,eta,mu)])

    return np.dot( Ns.T, dPHIs )


def Newton_iteration( xis, Ns, Xs, tol, max_it):
    "Newton iteration to compute the inverse transformation from \
    physical to master coordinates."
    
    R = Residual(xis[0], xis[1], xis[2], Ns, Xs)
    JF = -Jacobian(xis[0], xis[1], xis[2], Ns)
    U = -solve(JF,R)

    xis[0] = xis[0] + U[0][0]
    xis[1] = xis[1] + U[1][0]
    xis[2] = xis[2] + U[2][0]

    for i in range (max_it-1):
    
        R = Residual(xis[0], xis[1], xis[2], Ns, Xs)
        if norm(R) < tol:
            break

        JF = -Jacobian(xis[0], xis[1], xis[2], Ns)
        U = -solve(JF,R)

        xis[0] = xis[0] + U[0][0]
        xis[1] = xis[1] + U[1][0]
        xis[2] = xis[2] + U[2][0]
        

    return xis

def check_master_coord( xis ):
    "check if the master coordinates are in the domain [-1,1]^3"
    ## Caution!! some tolerane is needed for the fp arithmetic error

    if (xis[0]<= 1.0001) & (xis[0]>=-1.0001) &\
       (xis[1]<= 1.0001) & (xis[1]>=-1.0001) &\
       (xis[2]<= 1.0001) & (xis[2]>=-1.0001):
        return True
    else:
        return False




## Space integration of the error over the different meshes.
def space_integrate_error_2( xyzp_ref, xyzp_coa, nodes_ref, nodes_coa, \
                             elems_ref, elems_coa, max_depth, \
                             xyz_range ):
    "integrate the error between xyzp_ref (reference) and xyzp_coa (coarse)."
    # xyzp_ref: xyzp data on mesh 1
    # xyzp_coa: xyzp data on mesh 2
    # nodes_ref: nodes on mesh 1
    # nodes_coa: nodes on mesh 2
    # elems_ref: elements on mesh 1
    # elems_coa: elements on mesh 2

    xlen = xyz_range[1] - xyz_range[0]
    ylen = xyz_range[3] - xyz_range[2]
    zlen = xyz_range[5] - xyz_range[4]

    ## build octree of the corse mesh
    tree = []
    for i in range(len(elems_coa)):
        tree = build_octree(i, tree, elems_coa, nodes_coa, \
                            xyz_range, \
                            xlen, ylen, zlen, max_depth)
    
    Ex_error = 0.0
    Ey_error = 0.0
    Ez_error = 0.0
    Ep_error = 0.0
    
    # integrate over each element
    pr_c = 0
    for i in range(len(elems_ref)):
        # show progress
        if i > len(elems_ref)/10.0*pr_c:
            pr_c = pr_c+1
            sys.stdout.write('*')
            sys.stdout.flush()

        # pre-allocate the spaces
        xv = 8*[None]
        yv = 8*[None]
        zv = 8*[None]
        pv = 8*[None]

        for en in range(8):

            p = [nodes_ref[elems_ref[i][en]][0], \
                 nodes_ref[elems_ref[i][en]][1], \
                 nodes_ref[elems_ref[i][en]][2]]
            
            elem_n, xis = \
                inverse_physical_to_master( p, xyz_range, tree, \
                                            nodes_coa, elems_coa, \
                                            max_depth )

            xv[en], yv[en], zv[en], pv[en] = \
                get_xyzp_values( nodes_coa, elems_coa, xyzp_coa, elem_n, xis )


        # differences in xyzp values at nodes
        Ex = [xyzp_ref[elems_ref[i][0]][0] - xv[0], \
              xyzp_ref[elems_ref[i][1]][0] - xv[1], \
              xyzp_ref[elems_ref[i][2]][0] - xv[2], \
              xyzp_ref[elems_ref[i][3]][0] - xv[3], \
              xyzp_ref[elems_ref[i][4]][0] - xv[4], \
              xyzp_ref[elems_ref[i][5]][0] - xv[5], \
              xyzp_ref[elems_ref[i][6]][0] - xv[6], \
              xyzp_ref[elems_ref[i][7]][0] - xv[7]]
        Ey = [xyzp_ref[elems_ref[i][0]][1] - yv[0], \
              xyzp_ref[elems_ref[i][1]][1] - yv[1], \
              xyzp_ref[elems_ref[i][2]][1] - yv[2], \
              xyzp_ref[elems_ref[i][3]][1] - yv[3], \
              xyzp_ref[elems_ref[i][4]][1] - yv[4], \
              xyzp_ref[elems_ref[i][5]][1] - yv[5], \
              xyzp_ref[elems_ref[i][6]][1] - yv[6], \
              xyzp_ref[elems_ref[i][7]][1] - yv[7] ]
        Ez = [xyzp_ref[elems_ref[i][0]][2] - zv[0], \
              xyzp_ref[elems_ref[i][1]][2] - zv[1], \
              xyzp_ref[elems_ref[i][2]][2] - zv[2], \
              xyzp_ref[elems_ref[i][3]][2] - zv[3], \
              xyzp_ref[elems_ref[i][4]][2] - zv[4], \
              xyzp_ref[elems_ref[i][5]][2] - zv[5], \
              xyzp_ref[elems_ref[i][6]][2] - zv[6], \
              xyzp_ref[elems_ref[i][7]][2] - zv[7]] 
        Ep = [xyzp_ref[elems_ref[i][0]][3] - pv[0], \
              xyzp_ref[elems_ref[i][1]][3] - pv[1], \
              xyzp_ref[elems_ref[i][2]][3] - pv[2], \
              xyzp_ref[elems_ref[i][3]][3] - pv[3], \
              xyzp_ref[elems_ref[i][4]][3] - pv[4], \
              xyzp_ref[elems_ref[i][5]][3] - pv[5], \
              xyzp_ref[elems_ref[i][6]][3] - pv[6], \
              xyzp_ref[elems_ref[i][7]][3] - pv[7]] 


        Es_erros = integrate_hex8_xyzp( [nodes_ref[elems_ref[i][0]], \
                                         nodes_ref[elems_ref[i][1]], \
                                         nodes_ref[elems_ref[i][2]], \
                                         nodes_ref[elems_ref[i][3]], \
                                         nodes_ref[elems_ref[i][4]], \
                                         nodes_ref[elems_ref[i][5]], \
                                         nodes_ref[elems_ref[i][6]], \
                                         nodes_ref[elems_ref[i][7]]], \
                                        [Ex, Ey, Ez, Ep] )

        Ex_error = Ex_error + Es_erros[0]
        Ey_error = Ey_error + Es_erros[1]
        Ez_error = Ez_error + Es_erros[2]
        Ep_error = Ep_error + Es_erros[3]


    Ex_error = sqrt(Ex_error)
    Ey_error = sqrt(Ey_error)
    Ez_error = sqrt(Ez_error)
    Ep_error = sqrt(Ep_error)

    return (Ex_error, Ey_error, Ez_error, Ep_error)


def inverse_physical_to_master ( P, xyz_range, tree, nodes, elems, \
                                 max_depth ):
    "inverse transformation from phsyical coordinates to  element number \
    and master coordinates"    
    # P: physical coordinates
    # xyz_range: [xmin, xmax, ymin, ymax, zmin, zmax]

    xlen = xyz_range[1] - xyz_range[0]
    ylen = xyz_range[3] - xyz_range[2]
    zlen = xyz_range[5] - xyz_range[4]

    # get the traversal path to the branch 
    branch_chain = search_branch( P[0], P[1], P[2],  \
                                  xyz_range, \
                                  xlen, ylen, zlen, max_depth)

    # get candidate elements
    branch = tree
    if max_depth == 0:
        aux_max_depth = 1
    else:
        aux_max_depth = max_depth
    for dp in range(aux_max_depth):
        branch = branch[branch_chain[-1-dp]]

    if branch == []:
        print "branch empty, something went wrong!"
        print branch_chain
        print tree[0][0][4][0][0]

    # first guess for the master coordinates
    xis_guess = np.array([0.0,0.0,0.0])
    
    # allocate spaces
    Ns = np.empty([8,3])
    xv = np.empty([8,1])
    yv = np.empty([8,1])
    zv = np.empty([8,1])

    # tolerance and max iteration for the Newton iteration
    tol = 0.001
    max_it = 10

    for cn in range(len(branch)):
        for nn in range(8):
            Ns[nn] = ([nodes[elems[branch[cn]][nn]][0], \
                       nodes[elems[branch[cn]][nn]][1], \
                       nodes[elems[branch[cn]][nn]][2]])
        xis = Newton_iteration ( xis_guess, Ns, P, tol, max_it )

        # check if the master coordinates is in the legit domain
        if check_master_coord(xis):
            
            # return the element number and master coordinates
            return branch[cn], xis
            
    print "Cannot find the inverse in range!!!", P
    print "branch is:", branch

    for cn in range(len(branch)):
        for nn in range(8):
            Ns[nn] = ([nodes[elems[branch[cn]][nn]][0], \
                       nodes[elems[branch[cn]][nn]][1], \
                       nodes[elems[branch[cn]][nn]][2]])
        xis = Newton_iteration ( xis_guess, Ns, P, tol, max_it )
        print xis


def get_xyzp_values( nodes, elems, xyzp, elem_n, xis ):
    "get the values of x,y,z,p from the element number and \
    master coordinates."

    xv = []
    yv = []
    zv = []
    pv = []

    for en in range(8):
        # xv.append( nodes[elems[elem_n][en]][0] )
        # yv.append( nodes[elems[elem_n][en]][1] )
        # zv.append( nodes[elems[elem_n][en]][2] )

        # print elem_n

        xv.append (xyzp[elems[elem_n][en]][0])
        yv.append (xyzp[elems[elem_n][en]][1])
        zv.append (xyzp[elems[elem_n][en]][2])
        pv.append (xyzp[elems[elem_n][en]][3])

    px = interpolate(xis, xv)
    py = interpolate(xis, yv)
    pz = interpolate(xis, zv)
    pp = interpolate(xis, pv)

    return px, py, pz, pp


