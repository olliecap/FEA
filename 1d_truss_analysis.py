import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import inv

# element class makes the 1d beam and takes all necessary information. also makes local stiffness matrix to be used later

class element:
    def __init__(self, startNode, endNode, material, area, length):
        self.startNode = startNode
        self.endNode = endNode
        self.material = material
        self.area = area
        self.length = length
    def local_stiffness_matrix(self):
        k = (self.material * self.area) / self.length
        mat = k* np.array([[1, -1],
        [-1, 1]])
        return mat

# node class makes node for each point necessary. force and displacement depend on the boundary condition

class node:
    def __init__(self, coordinate, force, displacement):
        self.coordinate = coordinate
        self.force = force
        self.displacement = displacement
    def __str__(self):
        return f"Node {self.coordinate}\nForce: {self.force}\nDisplacement: {self.displacement}\n"

# makes a list with all the nodes using the node class

def nodeMaker():
    n = int(input("How many nodes are in your truss? "))
    print("Describe each node, from left to right.")
    nodes = []
    for i in range(n):
        force = input("What is the force acting on this node? Enter 'u' if unknown. ")
        displacement = input("What is the displacement on this node? Enter 'u' if unknown. ")
        print("\nNext node:\n")
        if force.lower() in ['u', 'unknown']:
            force_value = np.nan
        else:
            force_value = float(force)
        if displacement.lower() in ['u', 'unknown']:
            displacement_value = np.nan
        else:
            displacement_value = float(displacement)
        nodes.append(node(i, force_value, displacement_value))
    return nodes

# creates the elements and puts them all in a list, makes local stiffness matrix list as well

def element_maker(nodeList):
    print("Enter element values from left to right.\n")
    element_list = []
    for i in range(len(nodeList)-1):
        material = int(input("What is the E value of this element? "))
        area = int(input("What is the cross-sectional area of this element? "))
        length = int(input("What is the length of this element? "))
        element_list.append(element(i, i + 1, float(material), float(area), float(length)))
    LSM_list = []
    for s in element_list:
        LSM_list.append(s.local_stiffness_matrix())
    return element_list, LSM_list

# creates displacement matrix for solver

def displacement_matrix(nodeList):
    return [[node.displacement] for node in nodeList]

# creates force matrix for solver

def force_matrix(nodeList):
    return [[node.force] for node in nodeList]

# combines all the local matrices into one global stiffness matrix

def global_matrix(nodeList, matrix_list):
    glob = np.zeros((len(nodeList), len(nodeList)))
    for i, local in enumerate(matrix_list):
        start_node = i
        end_node = i + 1
        glob[start_node:start_node + 2, start_node:start_node + 2] += local
    return glob

# simplifies the gsm by inserting boundary conditions at fixed points

def boundary_conditions(nodeList, gsm, fMatrix):
    for i, node in enumerate(nodeList):
        if node.displacement == 0:
            for c, column in enumerate(gsm[i]):
                gsm[i][c] = 0
            for r, row in enumerate(gsm):
                gsm[r][i] = 0
            gsm[i][i] = 1
            fMatrix[i] = 0
    return gsm, fMatrix

"""For now, only solves the reaction forces in each end of the one element beam as of right now.
Takes the truss stiffness matrix, displacements from both fixed ends, and solve for the reaction
forces. Displacement and force are taken from nodeList, k has local_stiffness_matrix """

"""
MAKE A FOR LOOP THAT PARTITIONS BOTH THE STIFFNESS MATRIX WITH BOUNDARY conditions
AND THE STIFFNESS MATRIX WITHOUT BOUNDARY CONDITIONS, COMMENT EVERYWHERE AND FIX NAMES AS WELL
"""
# bck = boundary conditioned k, ogk = orignal k matrix
def solver(f, bck, ogk, u):
    u_displacement = []
    k_displacement = []
    partition_k_displacement = []
    partition_k_force = []
    for i, displacement in enumerate(u):
        if np.isnan(displacement):
            u_displacement.append(i)
            partition_k_force.append(f[i])
        else:
            k_displacement.append(i)
            partition_k_displacement.append(displacement)
    partition_k_force = np.array(partition_k_force)
    partition_k_displacement = np.array(partition_k_displacement)

    # partition the gsm into 4 parts, known known, unknown known, known unknown, and unknown unknown, based on what displacements we know
    k_uu = k[np.ix_(u_displacement, u_displacement)]
    k_uk = k[np.ix_(u_displacement, k_displacement)]
    k_ku = k[np.ix_(k_displacement, u_displacement)]
    k_kk = k[np.ix_(k_displacement, k_displacement)]

    # solves for unknown displacements
    partition_u_displacement = spsolve(k_uu, partition_k_force - k_uk @ partition_k_displacement)

    # solves for unknown forces
    partition_u_force = k_ku @ partition_u_displacement + k_kk @ partition_k_displacement

    # adds the unknown forces and displacments back to the orignal matrices
    z = 0
    for i, displacement in enumerate(u):
        if np.isnan(displacement):
            u[i] = partition_u_displacement[z]
            z += 1
    z = 0
    for i, force in enumerate(f):
        if np.isnan(force):
            u[z] = partition_u_force[z]
            z += 1
    return f , u

def main():
    allTheNodes = nodeMaker()
    elements, lsms = element_maker(allTheNodes)
    umat = displacement_matrix(allTheNodes)
    fmat = force_matrix(allTheNodes)
    gm = global_matrix(allTheNodes, lsms)
    gsm, fmat = boundary_conditions(allTheNodes, gm, fmat)
    print(gsm)
    solver(fmat, gsm, gm, umat)

if __name__ == '__main__':
    main()
