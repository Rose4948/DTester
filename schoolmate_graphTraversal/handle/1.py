def LCS_length_3(x, y):
    xlen = len(x);
    ylen = len(y);
    tmp_1 = [[0 for i in range(ylen + 1)] for j in range(2)]
    tmp_2 = [[0 for i in range(ylen + 1)] for j in range(xlen + 1)]
    for i in range(1, xlen + 1):
        for j in range(1, ylen + 1):
            if (x[i - 1] == y[j - 1]):
                tmp_1[i % 2][j] = tmp_1[(i - 1) % 2][j - 1] + 1;
                tmp_2[i][j] = 0
            else:
                if (tmp_1[i % 2][j - 1] >= tmp_1[(i - 1) % 2][j]):
                    tmp_1[i % 2][j] = tmp_1[i % 2][j - 1]
                    tmp_2[i][j] = 1
                else:
                    tmp_1[i % 2][j] = tmp_1[(i - 1) % 2][j]
                    tmp_2[i][j] = -1
    return tmp_1, tmp_2
def LCS_Len_And_Index(x, y, tmp_2):
    i = len(x)
    j = len(y)
    index_x=i/2
    index_y=j/2
    k=0
    while (i > 0 and j > 0):
        if (tmp_2[i][j] == 0):
            index_x=i-1
            index_y=j-1
            k+=1
            i -= 1
            j -= 1
        elif tmp_2[i][j] == 1:
            j -= 1
        elif tmp_2[i][j] == -1:
            i -= 1
    return k,index_x,index_y

#there are LCS Sequce and index of x,y
def LCS_print(x, y, tmp_2):
    result = []
    i = len(x)
    j = len(y)
    k = 0
    index_x=i/2
    index_y=j/2
    while (i > 0 and j > 0):
        if (tmp_2[i][j] == 0):
            result.append(x[i - 1])
            index_x=i-1
            index_y=j-1
            k = k + 1
            i = i - 1
            j = j - 1
        elif tmp_2[i][j] == 1:
            j = j - 1
        elif tmp_2[i][j] == -1:
            i = i - 1
    return result,index_x,index_y

def test_3(x, y):
    max_len, tmp_2 = LCS_length_3(x, y)
    LCS_Len,index_x,index_y = LCS_Len_And_Index(x, y, tmp_2)
    return LCS_Len,index_x,index_y

def main():

    x=['T12', 'T14', 'T23','T23', 'T15', 'T19','T15']
    y=['T10', 'T11', 'T14', 'T15', 'T16', 'T19']




if __name__ == '__main__':  # not execute when import as a module
    main()
if __name__ == '__main__':  # not execute when import as a module
    main()