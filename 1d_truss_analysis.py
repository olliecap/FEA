import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

"""
My very first FEM project :). Takes a 1d beam with theoretically infinite nodes and spits out the displacement and forces at each node. Very simple project but a
good starting one to learn both FEA and Python, as someone who has limited knowledge in both subjects
"""

# element class makes the 1d beam and takes all necessary information. also makes local stiffness matrix to be used later
class element:
    def __init__(self, startNode, endNode, material, area, length):
        self.startNode = startNode
        self.endNode = endNode
        self.material = material
        self.area = area
        self.length = length

    # local stiffness matrix (lsm) for each node. will later be added to a global stiffness matrix for solver
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

    # how many nodes?
    n = int(input("How many nodes are in your truss? "))

    # based on node count, iterates through number and creates node using class. appends it to nodes list
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

# creates the elements and puts them all in a list, makes local stiffness matrix list as well. based off inputs
def element_maker(nodeList):

    # takes information for each element. number of elements is n-1 nodes
    print("Enter element values from left to right.\n")
    element_list = []
    for i in range(len(nodeList)-1):
        material = int(input("\nWhat is the E value of this element? "))
        area = int(input("\nWhat is the cross-sectional area of this element? "))
        length = int(input("\nWhat is the length of this element? "))
        element_list.append(element(i, i + 1, float(material), float(area), float(length)))

    # put all the elements in a list
    LSM_list = []
    for s in element_list:
        LSM_list.append(s.local_stiffness_matrix())
    return element_list, LSM_list

# creates displacement vector for solver
def displacement_matrix(nodeList):
    return [[node.displacement] for node in nodeList]

# creates force vector for solver
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
def boundary_conditions(nodeList, gm, fMatrix):
    gsm = gm.copy()
    for i, node in enumerate(nodeList):
        if node.displacement == 0:
            for c, column in enumerate(gsm[i]):
                gsm[i][c] = 0
            for r, row in enumerate(gsm):
                gsm[r][i] = 0
            gsm[i][i] = 1
            fMatrix[i] = 0
    return gsm, fMatrix

# bck = boundary conditioned k, ogk = orignal k matrix
def solver(f, bck, ogk, u):

    # partitions the u and k vectors. splits them into known and unknown submatrices
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

    # turns the lists into arrays for matrix multiplication
    partition_k_force = np.array(partition_k_force)
    partition_k_displacement = np.array(partition_k_displacement)

    # partition bck into 4 parts, known known, unknown known, known unknown, and unknown unknown, based on what displacements we know
    k_uu = csr_matrix(bck[np.ix_(u_displacement, u_displacement)])
    k_uk = csr_matrix(bck[np.ix_(u_displacement, k_displacement)])
    k_ku = csr_matrix(bck[np.ix_(k_displacement, u_displacement)])
    k_kk = csr_matrix(bck[np.ix_(k_displacement, k_displacement)])

    # solves for unknown displacements
    partition_u_displacement = spsolve(k_uu, partition_k_force - k_uk @ partition_k_displacement)

    # partition ogk
    k_uu = ogk[np.ix_(u_displacement, u_displacement)]
    k_uk = ogk[np.ix_(u_displacement, k_displacement)]
    k_ku = ogk[np.ix_(k_displacement, u_displacement)]
    k_kk = ogk[np.ix_(k_displacement, k_displacement)]

    # turns both partitions into 1x2 matrices, so solving for the forces always gives a 1x2 matrix
    partition_k_displacement = partition_k_displacement.flatten()
    partition_u_displacement = partition_u_displacement.flatten()

    # solves for unknown forces
    partition_u_force = k_ku @ partition_u_displacement + k_kk @ partition_k_displacement

    # adds the unknown forces and displacments back to the orignal matrices
    t = 0
    z = 0
    for i, displacement in enumerate(u):
        if np.isnan(displacement):
            try:
                u[i] = partition_u_displacement[t]
            except:
                pass
            t = 1
        else:
            try:
                f[i] = partition_u_force[z]
            except:
                pass
            z += 1

    # there are some 1 element lists in the f and u lists, this turns them into floats
    for i, val in enumerate(u):
        if isinstance(val, list):
            u[i] = val[0]

    for i, val in enumerate(f):
        if isinstance(val, list):
            f[i] = val[0]

    # printing the outputs
    print("Displacement:")
    for i, val in enumerate(u):
        print(f"Node {i}: {val}")

    print("Force")
    for i, val in enumerate(f):
        print(f"Node {i}: {val}")

    return f , u

def main():
    allTheNodes = nodeMaker()
    elements, lsms = element_maker(allTheNodes)
    umat = displacement_matrix(allTheNodes)
    fmat = force_matrix(allTheNodes)
    gm = global_matrix(allTheNodes, lsms)
    gsm, fmat = boundary_conditions(allTheNodes, gm, fmat)
    solver(fmat, gsm, gm, umat)

if __name__ == '__main__':
    main()
