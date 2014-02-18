import numpy as np 

def read_nodes_elems( filename ):

    ## get number of nodes and elements
    n_flag = 0
    e_flag = 0
    n_nodes = 0
    n_elems = 0
    with open(filename) as f:
        for line in f:
            if(line.find("</Nodes>")!=-1):
                n_flag = 0
            elif line.find("</Elements>") != -1:
                e_flag = 0

            if n_flag==1:
                n_nodes = n_nodes + 1
            if e_flag==1:
                n_elems = n_elems + 1

            if(line.find("<Nodes>")!=-1):
                n_flag = 1
            if(line.find("<Elements>")!=-1):
                e_flag = 1

    ### LOAD NODE POSITIONS ###
    f = open(filename)
    content = f.readlines()
    f.close()
    len_file = len(content)
            

    nodes = np.empty([n_nodes, 3], dtype=float)
    elems = np.empty([n_elems,8], dtype=int)

    n_flag = 0
    e_flag = 0
    for i in range(len_file):
        if content[i].find("</Nodes>") != -1:
            n_flag = 0
        elif content[i].find("</Elements>") != -1:
            e_flag = 0

        ## Add Nodes
        if n_flag > 0:
            node_sp = content[i].split('>')[1].split('<')[0].split(',')
            nodes[n_flag-1] = [float(node_sp[0]), float(node_sp[1]), float(node_sp[2])]
            n_flag = n_flag+1

        ## Add Elements
        if e_flag > 0:
            elem_sp = content[i].split('>')[1].split('<')[0].split(',')
            ## "-1" for these because the index starts from 0
            elems[e_flag-1] = [int(elem_sp[0])-1, int(elem_sp[1])-1, 
                    int(elem_sp[2])-1, int(elem_sp[3])-1,
                    int(elem_sp[4])-1, int(elem_sp[5])-1,
                    int(elem_sp[6])-1, int(elem_sp[7])-1]
            e_flag = e_flag+1

        if content[i].find("<Nodes>") != -1:
            n_flag = 1
        elif content[i].find("<Elements>") != -1:
            e_flag = 1

    return (nodes, elems)


def read_xyzp( filename, n_nodes, end_step ):
    "read the x, y, z, p values from the result file"
    # end_step: ending step of the simulation
    # end_step=102 for mesh erorr analysis
    # end_step=-1 for time error analysis
    
    f = open(filename)
    content = f.readlines()
    f.close()
    len_file = len(content)

    xyzps_set = []

    # first get the number of nodes
    for i in range(len_file):
        if content[i].find('*Step  = 2') != -1:
            n_nodes = int(content[i-1].split(' ')[0])
                        
    xyzps = np.empty([n_nodes, 4], dtype=float)

    for i in range(len_file):            
        if ((content[i].find('*Time') != -1) & \
            (content[i-1].find(str(end_step)) != -1)) | \
            ((content[i].find('*Time  = 2') != -1) & (end_step==-1)) :

            # for the time step number 'Step', read the x,y,z,p
            Step = int(content[i-1].split('=')[1])
            # xyzps = [[0.0 for ii in range(4)] for jj in range(n_nodes)]
            for k in range(n_nodes):
                xyzp_split = content[i+2+k].split(' ')
                xyzps[k][0] = float(xyzp_split[1])
                xyzps[k][1] = float(xyzp_split[2])
                xyzps[k][2] = float(xyzp_split[3])
                xyzps[k][3] = float(xyzp_split[4])

            xyzps_set.append (xyzps)

    return xyzps_set

def write_geo_file( sim_number, layers ):

    f1 = open("box_tmpl.geo", 'r')
    f2 = open('box.geo', 'w')
    for line in f1:
        f2.write(line.replace('layer = 16', 'layer = %i' %layers))
    f1.close()
    f2.close()
