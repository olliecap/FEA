# Matrix Operations Library
def get_matrix(row, col):
    matrix = []
    print("Enter matrix values (separated by space, rows by new line)")
    for i in range(row):
        rowString = input()
        rowStringList = rowString.split()
        while len(rowStringList) != col:
            print("Invalid num of columns, columns must be " + str(col))
            rowString = input()
            rowStringList = rowString.split()
        r = []
        for num in rowStringList:
            r.append(int(num))
        matrix.append(r)
    return matrix

def printMat(matrix):
    for r in range(len(matrix)):
        for c in range(len(matrix[0])):
            print(str(matrix[r][c]) + " ", end='')
        print()

def dim(number):
    bothMatrix = []
    for i in range(number):
        matDim = input("What are the dimensions of your matrix (ex: 2x3) ")
        dimensions = matDim.split("x")
        row = int(dimensions[0].strip())
        col = int(dimensions[1].strip())
        first = get_matrix(row, col)
        bothMatrix.append(first)
    return bothMatrix
# multiply function, mostly to boost my python knowledge
def multiply():
    result = []
    first, second = dim(2)
    check = False
    while check == False:
        if len(first[0]) != len(second):
            print("Number of rows in first matrix must equal number of columns in second matrix.\n")
            first, second = dim()
        else:
            check = True
    try:
        for r in range(len(first)):
            row = []
            for t in range(2):
                integer = []
                for c in range(len(second)):
                    integer.append(first[r][c]*second[c][t])
                row.append(sum(integer))
            result.append(row)
        print(result)
    except:
        pass
    # print result so it looks like a matrix
    for t in range(len(result)):
        if t%len(first) == 0 and t != 0:
            print()
        print(str(result[t]) + " ", end='')
    return result
# adding function, mostly to boost my python knowledge
def add():
    first, second = dim(2)
    result = []
    for r in range(len(first)):
        resultRow = []
        for c in range(len(first[0])):
            resultRow.append(first[r][c] + second[r][c])
        result.append(resultRow)
    return printMat(result), result
# transpose function, mostly to boost my python knowledge
def transpose(mat):
    newMat = []
    for c in range(len(mat[0][0])):
        row = []
        for r in range(len(mat[0])):
            row.append(mat[0][r][c])
        newMat.append(row)
    printMat(newMat)
    return newMat

"""Gauss-Seidel method. sys is a list. inside is more lists, each being its own
equation. The first value is the coefficient of x1, second is x2 coefficient, etc.
The last one is the constant for the equation. could be useful for solving long,
diagonally domininant martices"""
def seidel(sys, tolerance = 0.001, max_iterations = 1000):
    # Guessed values
    values = [0 for eq in range(len(sys[0])-1)]
    # Start of Guass-Seidel method
    for iteration in range(max_iterations):
        # oldVals will be used to check if method has converged
        oldVals = values.copy()
        for i in range(len(values)):
            # take the constant out of the eq so it doesn't get messed with the coefficients
            cons = sys[i].pop()
            newCoeffs = 0
            for c in range(len(sys[i])):
                # check if we hit the diagonal element, aka the x that we are solving for. if not, calculate the new x value to be added to the constant
                if c != i:
                    newCoeffs += -sys[i][c] * values[c]
            # change value of the x's, which will be used in the next iteration
            values[i] = (cons + newCoeffs) / sys[i][i]
            sys[i].append(cons)
        # if we are within specified tolerance, quit the function
        if all(abs(values[i] - oldVals[i]) <= tolerance for i in range(len(sys))):
            print(values)
            return values
        print(values)
    # didn't reach convergence within specified number of iterations
    print("Didn't converge")

def main():
    validInput = True
    while validInput:
        cmd = input("Enter + to add, * to multiply, t to transpose, s for seidel: ")
        if cmd == "+":
            validInput = False
            add()
        elif cmd == "*":
            validInput = False
            multiply()
        elif cmd == "t":
            validInput = False
            transpose(dim(1))
        elif cmd == "s":
            validInput = False
            seidel([[4, 1, 1, 6], [2, 6, 1, 9], [1, 1, 5, 7]])
        else:
            print("Invalid input\n")
if __name__ == "__main__":
    main()
