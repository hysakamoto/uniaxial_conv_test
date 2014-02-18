from subprocess import call
import xml.etree.ElementTree as etree
# from lxml import etree

##### simulation
def simulation(id, sim_name, name, pressure_top, pressure_side, \
               time_steps, min_dtmax, step_size, \
               phi0, density, \
               c1, c2, k, \
               ksi, beta, \
               g1, t1, \
               perm, \
               sliding_penalty\
):

    f = open('box.msh')
    content = f.readlines()
    f.close()

    # output
    call(["mkdir", "-p", sim_name])
    feb_file = '{0:s}/{1:s}.feb'.format(sim_name,name)
    f = open(feb_file, 'w')

    # scale
    # sc = 1e-6
    sc = 1

    len_file = len(content)
    node_start = len_file
    node_end = -1
    elem_start = len_file
    elem_end = -1

    Nodes = []

    hex8_cell = []

    quad4_celltop =[]
    quad4_cellbottom = []
    quad4_xz = []
    quad4_yz = []
    quad4_cellside = []

    node_topplane = []
    node_xyplane = []
    node_xzplane = []
    node_yzplane = []
    node_side = []

    for i in range(len_file):
        # get number of nodes
        # get the start and end line numbers of node data
        if content[i].find("$Nodes") != -1:
            n_nodes = int(content[i+1])
            node_start = i+2
            node_end = i+2+n_nodes

        # collect node information
        if (i>=node_start) & (i<node_end):
            node_list = content[i].split('\n')[0].split(' ')
            node_info = []
            node_info.append(int(node_list[0]))
            node_info.append(float(node_list[1])*sc)
            node_info.append(float(node_list[2])*sc)
            node_info.append(float(node_list[3])*sc)
            Nodes.append(node_info)

        # get number of elements
        # get the start and end line numbers of element data
        if content[i].find("$Elements") != -1:
            n_elem = int(content[i+1])
            elem_start = i+2
            elem_end = i+2+n_elem

        # collect element information
        if (i>=elem_start) & (i<elem_end):
            elem_list = content[i].split('\n')[0].split(' ')
            elem_info = []
            elem_info.append(int(elem_list[0]))
            elem_info.append(int(elem_list[1]))
            elem_info.append(int(elem_list[2]))
            elem_info.append(int(elem_list[3]))
            elem_info.append(int(elem_list[4]))

            # hex8 elements
            if (elem_info[1] == 5):
                elem_info.append(int(elem_list[5]))
                elem_info.append(int(elem_list[6]))
                elem_info.append(int(elem_list[7]))
                elem_info.append(int(elem_list[8]))
                elem_info.append(int(elem_list[9]))
                elem_info.append(int(elem_list[10]))
                elem_info.append(int(elem_list[11]))
                elem_info.append(int(elem_list[12]))
                # cell
                if  (elem_info[3]==1):
                    hex8_cell.append(elem_info)

            # quad4 elements
            if (elem_info[1] == 3):
                elem_info.append(int(elem_list[5]))
                elem_info.append(int(elem_list[6]))
                elem_info.append(int(elem_list[7]))
                elem_info.append(int(elem_list[8]))
                # cell bottom surface
                if (elem_info[3]==1):
                    quad4_cellbottom.append(elem_info)
                # cell top surface
                elif (elem_info[3]==2):
                    quad4_celltop.append(elem_info)
                elif (elem_info[3]==3):
                    quad4_xz.append(elem_info)
                elif (elem_info[3]==4):
                    quad4_yz.append(elem_info)
                # cell side
                elif (elem_info[3]==5) | (elem_info[3]==6):
                    quad4_cellside.append(elem_info)

    # collect and order the top
    for i in range(len(quad4_celltop)):
        node_topplane.append(quad4_celltop[i][5])
        node_topplane.append(quad4_celltop[i][6])
        node_topplane.append(quad4_celltop[i][7])
        node_topplane.append(quad4_celltop[i][8])
    # remove duplicates and sort
    node_topplane=list(set(node_topplane))
    node_topplane.sort()

    # collect and order the bottom
    for i in range(len(quad4_cellbottom)):
        node_xyplane.append(quad4_cellbottom[i][5])
        node_xyplane.append(quad4_cellbottom[i][6])
        node_xyplane.append(quad4_cellbottom[i][7])
        node_xyplane.append(quad4_cellbottom[i][8])
    # remove duplicates and sort
    node_xyplane=list(set(node_xyplane))
    node_xyplane.sort()

    # collect and order the bottom
    for i in range(len(quad4_xz)):
        node_xzplane.append(quad4_xz[i][5])
        node_xzplane.append(quad4_xz[i][6])
        node_xzplane.append(quad4_xz[i][7])
        node_xzplane.append(quad4_xz[i][8])
    # remove duplicates and sort
    node_xzplane=list(set(node_xzplane))
    node_xzplane.sort()

    # collect and order the bottom
    for i in range(len(quad4_yz)):
        node_yzplane.append(quad4_yz[i][5])
        node_yzplane.append(quad4_yz[i][6])
        node_yzplane.append(quad4_yz[i][7])
        node_yzplane.append(quad4_yz[i][8])
    # remove duplicates and sort
    node_yzplane=list(set(node_yzplane))
    node_yzplane.sort()

    # collect and order the side
    for i in range(len(quad4_cellside)):
        node_side.append(quad4_cellside[i][5])
        node_side.append(quad4_cellside[i][6])
        node_side.append(quad4_cellside[i][7])
        node_side.append(quad4_cellside[i][8])
    # remove duplicates and sort
    node_side=list(set(node_side))
    node_side.sort()


    #####################################################################
    ##                Now write the .feb file                          ##
    #####################################################################

    # header
    root = etree.Element("febio_spec", version="1.2")

    # Global
    Globals_xml = etree.SubElement(root, "Globals")

    ## Constants
    Constants_xml = etree.SubElement(Globals_xml, "Constants")
    etree.SubElement(Constants_xml, "T").text = "0"
    etree.SubElement(Constants_xml, "R").text = "0"
    etree.SubElement(Constants_xml, "Fc").text = "0"

    # Material
    Material_xml = etree.SubElement(root, "Material")

    material_xml = etree.SubElement(Material_xml, "material", id = "1", name = "cell", type = "biphasic")
    etree.SubElement(material_xml, "phi0").text = "{0:e}".format(phi0)
    solid_xml = etree.SubElement(material_xml, "solid", type="Mooney-Rivlin")
    etree.SubElement(solid_xml, "density").text = "{0:e}".format(1)
    etree.SubElement(solid_xml, "c1").text = "{0:e}".format(c1)
    etree.SubElement(solid_xml, "c2").text = "{0:e}".format(0)
    etree.SubElement(solid_xml, "k").text = "{0:e}".format(k)
    etree.SubElement(solid_xml, "laugon").text = "1"
    etree.SubElement(solid_xml, "atol").text = "0.05"
    permeability_xml = etree.SubElement(material_xml, "permeability", type="perm-const-iso")
    etree.SubElement(permeability_xml, "perm").text = "{0:e}".format(perm)

    ## cell
    # material_xml = etree.SubElement(Material_xml, "material", id = "1", name = "cell", type = "Mooney-Rivlin")

    # # etree.SubElement(material_xml, "density").text = "1"
    # etree.SubElement(material_xml, "c1").text = "{0:d}".format(c1)
    # etree.SubElement(material_xml, "c2").text = "{0:d}".format(0)
    # etree.SubElement(material_xml, "k").text = "{0:d}".format(k)
    # etree.SubElement(material_xml, "laugon").text = "1"
    # etree.SubElement(material_xml, "atol").text = "0.05"

    # Geometry
    Geometry_xml = etree.SubElement(root, "Geometry")

    ## Nodes
    Nodes_xml = etree.SubElement(Geometry_xml, "Nodes")
    for i in range(len(Nodes)):
        etree.SubElement(Nodes_xml, "node", id="{0:d}".format(Nodes[i][0])).text = "{0:9.7e}, {1:9.7e}, {2:9.7e}".format(Nodes[i][1], Nodes[i][2], Nodes[i][3])

    ## Elements
    Elements_xml = etree.SubElement(Geometry_xml, "Elements")
    elem_id = 1

    ### hex8
    for i in range(len(hex8_cell)):
        hex8_cell[i][0]=elem_id # fix element id
        elem_id+=1;
        hex8_cell[i][3]=1 # fix material id
        etree.SubElement(Elements_xml, "hex8",
                         id="{0:d}".format(hex8_cell[i][0]),
                         mat="{0:d}".format(hex8_cell[i][3])).text = "{0:5d}, {1:5d}, {2:5d}, {3:5d}, {4:5d}, {5:5d}, {6:5d}, {7:5d}".format(hex8_cell[i][5], hex8_cell[i][6], hex8_cell[i][7], hex8_cell[i][8], hex8_cell[i][9], hex8_cell[i][10], hex8_cell[i][11], hex8_cell[i][12])

    # Boundary
    Boundary_xml = etree.SubElement(root, "Boundary")

    ## fixed boundary on xy-plane (no z-displacement)
    fix_xml = etree.SubElement(Boundary_xml, "fix")
    for i in range(len(node_xyplane)):
        etree.SubElement(fix_xml, "node", id="{0:d}".format(node_xyplane[i]), bc="z")

    ## fixed boundary on xz-plane (no y-displacement)
    fix_xml = etree.SubElement(Boundary_xml, "fix")
    for i in range(len(node_xzplane)):
        etree.SubElement(fix_xml, "node", id="{0:d}".format(node_xzplane[i]), bc="y")

    ## fixed boundary on yz-plane (no x-displacement)
    fix_xml = etree.SubElement(Boundary_xml, "fix")
    for i in range(len(node_yzplane)):
        etree.SubElement(fix_xml, "node", id="{0:d}".format(node_yzplane[i]), bc="x")

    ## fixed pressure
    fix_xml = etree.SubElement(Boundary_xml, "fix")
    for i in range(len(node_side)):
        etree.SubElement(fix_xml, "node", id="{0:d}".format(node_side[i]), bc="p")


    # ## prescribed pressre on xy-plane
    # prescribe_xml = etree.SubElement(Boundary_xml, "prescribe")
    # for i in range(len(node_yzplane)):
    #     etree.SubElement(prescribe_xml, "node", id="{0:d}".format(node_yzplane[i]), bc="p", lc="3").text = "1"

    # ## prescribed pressre on top-plane
    # prescribe_xml = etree.SubElement(Boundary_xml, "prescribe")
    # for i in range(len(node_topplane)):
    #     etree.SubElement(prescribe_xml, "node", id="{0:d}".format(node_topplane[i]), bc="p", lc="3").text = "-1"


    # Load Data
    LoadData_xml = etree.SubElement(root, "LoadData")

    ### saving step
    loadcurve_xml = etree.SubElement(LoadData_xml, "loadcurve", id="1", type="linear", extend="constant")
    for i in range(100):
        etree.SubElement(loadcurve_xml, "loadpoint").text = "{0:e}, {1:e}".format(float(i)*0.1, min_dtmax)
        
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "0.1, 0"
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "0.2, 0"
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "0, 0"
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "0, 0"
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "0, 0"
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "0, 0"
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "0, 0"
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "1.0, {0:e}".format(min_dtmax)
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "2.0, {0:e}".format(min_dtmax)
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "3.0, {0:e}".format(min_dtmax)
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "4.0, {0:e}".format(min_dtmax)
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "5.0, {0:e}".format(min_dtmax)
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "6.0, {0:e}".format(min_dtmax)
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "7.0, {0:e}".format(min_dtmax)
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "8.0, {0:e}".format(min_dtmax)
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "9.0, {0:e}".format(min_dtmax)
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "10.0, {0:e}".format(min_dtmax)
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "11.0, {0:e}".format(min_dtmax)
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "12.0, {0:e}".format(min_dtmax)
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "13.0, {0:e}".format(min_dtmax*2)
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "14.0, {0:e}".format(min_dtmax*2)
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "15.0, {0:e}".format(min_dtmax*2)
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "16.0, {0:e}".format(min_dtmax*2)
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "17.0, {0:e}".format(min_dtmax*2)
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "18.0, {0:e}".format(min_dtmax*2)
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "19.0, {0:e}".format(min_dtmax*2)
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "20.0, {0:e}".format(min_dtmax*2)
    # etree.SubElement(loadcurve_xml, "loadpoint").text = "21.0, {0:e}".format(min_dtmax*2)

    ### traction application: top
    loadcurve_xml = etree.SubElement(LoadData_xml, "loadcurve", id="2", type="smooth")
    etree.SubElement(loadcurve_xml, "loadpoint").text = "0, 0"
    etree.SubElement(loadcurve_xml, "loadpoint").text = "1.0, {0:9.7e}".format(pressure_top)

    ### pressure application: top
    loadcurve_xml = etree.SubElement(LoadData_xml, "loadcurve", id="3", type="smooth")
    etree.SubElement(loadcurve_xml, "loadpoint").text = "0, 0"
    etree.SubElement(loadcurve_xml, "loadpoint").text = "1.0, {0:9.7e}".format(100.0)

    # Output
    Output_xml = etree.SubElement(root, "Output")
    plotfile_xml = etree.SubElement(Output_xml, "plotfile", type="febio")
    etree.SubElement(plotfile_xml, "var", type="displacement")
    etree.SubElement(plotfile_xml, "var", type="velocity")
    etree.SubElement(plotfile_xml, "var", type="effective fluid pressure")
    etree.SubElement(plotfile_xml, "var", type="fluid flux")
    etree.SubElement(plotfile_xml, "var", type="stress")
    etree.SubElement(plotfile_xml, "var", type="relative volume")

    logfile_xml = etree.SubElement(Output_xml, "logfile")
    # etree.SubElement(logfile_xml, "node_data", data="x;y;z;p", file="{0:s}.dat".format(name)).text="1:8:1"

    etree.SubElement(logfile_xml, "node_data", data="x;y;z;p", file="{0:s}/{1:s}.dat".format(sim_name,name))
    # etree.SubElement(logfile_xml, "rigid_body_data", data="Fz").text = "3"


    # Step01
    Step_xml = etree.SubElement(root, "Step", name="Step01")
    # etree.SubElement(Step_xml, "Module", type="biphasic")
    etree.SubElement(Step_xml, "Module", type="biphasic")
    Control_xml = etree.SubElement(Step_xml, "Control")

    etree.SubElement(Control_xml, "time_steps").text = "{0:d}".format(time_steps)
    etree.SubElement(Control_xml, "step_size").text = "%e"%step_size
    etree.SubElement(Control_xml, "max_refs").text = "15"
    etree.SubElement(Control_xml, "max_ups").text = "10"
    etree.SubElement(Control_xml, "dtol").text = "0.001"
    etree.SubElement(Control_xml, "etol").text = "0.01"
    etree.SubElement(Control_xml, "rtol").text = "0"
    # etree.SubElement(Control_xml, "ptol").text = "0.01"
    etree.SubElement(Control_xml, "lstol").text = "0.9"

    # time_stepper_xml = etree.SubElement(Control_xml, "time_stepper")
    # etree.SubElement(time_stepper_xml, "dtmin").text = "0.5"
    # etree.SubElement(time_stepper_xml, "dtmax").text = "{0:e}".format(min_dtmax)
    # etree.SubElement(time_stepper_xml, "dtmax", lc="1").text = "{0:e}".format(min_dtmax)
    # etree.SubElement(time_stepper_xml, "max_retries").text = "10"
    # etree.SubElement(time_stepper_xml, "opt_iter").text = "10"

    etree.SubElement(Control_xml, "plot_level").text = "PLOT_MUST_POINTS"
    etree.SubElement(Control_xml, "print_level").text = "PRINT_MAJOR_ITRS"

    ## Loads
    Loads_xml = etree.SubElement(Step_xml, "Loads")

    ### traction: top
    pressure_xml = etree.SubElement(Loads_xml, "traction")
    for i in range(len(quad4_celltop)):
        etree.SubElement(pressure_xml, "quad4", id="{0:d}".format(i+1), lc="2", tx="0", ty="0", tz="1").text = "{0:d}, {1:d}, {2:d}, {3:d}".format(quad4_celltop[i][5], quad4_celltop[i][6], quad4_celltop[i][7], quad4_celltop[i][8])

    ### fluidflux: symmetry planes
    fluidflux_xml = etree.SubElement(Loads_xml, "fluidflux", type="nonlinear", flux="fluid")
    for i in range(len(quad4_cellbottom)):
        etree.SubElement(fluidflux_xml, "quad4", id="{0:d}".format(i+1), lc="2", scale="0").text = "{0:d}, {1:d}, {2:d}, {3:d}".format(quad4_cellbottom[i][5], quad4_cellbottom[i][6], quad4_cellbottom[i][7], quad4_cellbottom[i][8])
    for i in range(len(quad4_xz)):
        etree.SubElement(fluidflux_xml, "quad4", id="{0:d}".format(i+1), lc="2", scale="0").text = "{0:d}, {1:d}, {2:d}, {3:d}".format(quad4_xz[i][5], quad4_xz[i][6], quad4_xz[i][7], quad4_xz[i][8])
    for i in range(len(quad4_yz)):
        etree.SubElement(fluidflux_xml, "quad4", id="{0:d}".format(i+1), lc="2", scale="0").text = "{0:d}, {1:d}, {2:d}, {3:d}".format(quad4_yz[i][5], quad4_yz[i][6], quad4_yz[i][7], quad4_yz[i][8])

    # get the string and write to file
    xml_str = etree.tostring(root, encoding='ISO-8859-1')
    # xml_str = etree.tostring(root)
    # f.write(xml_str)
    indent(root)

    # wrap the root inside ElementTree instance
    # tree = etree.ElementTree(root)
    # tree.write('output.xml', encoding='ISO-8859-1', xml_declaration='True')

    # root_parsed = etree.parse(feb_file).getroot()
    # indent(root_parsed)

    # etree.dump(root_parsed)

    # xml declaration
    f.write("""<?xml version="1.0" encoding="ISO-8859-1"?>
""")
    f.write(etree.tostring(root))

    # close the file
    f.close()

def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i




