import matplotlib.pyplot as plt

def main():
    M = list()
    # with open(r"/Users/h1007185/workspace/Rust/jm_lib_par/res/py/matrix.txt") as f:
    with open(r"res/py/matrix2.txt") as f:
        lines = f.readlines()
        for line in lines:
            # print(line.strip())
            row = list()
            for s in line.strip().split(" "):
                row.append(float(s))
            M.append(row)

    plt.matshow(M, interpolation=None, cmap="gray", vmin=0, vmax=1)
    plt.show()

    ## Symmetric check
    # m = len(M)
    # for i in range(0, m):
    #     for j in range(0, m):
    #         if M[i][j] != M[j][i]:
    #             print("{}, {}".format(i, j))

if __name__ == "__main__":
    main()