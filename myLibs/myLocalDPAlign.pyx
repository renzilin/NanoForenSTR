from collections import namedtuple

global F
F = namedtuple('F', ('score', 'pointer'))

def init_array(x, y):
    array = [[0] * (y) for _ in range(x)]
    array[0][0] = F(0, None)
    for j in range(1, y):
        array[0][j] = F((-5)*j, [0, j-1])
    for i in range(1, x):
        array[i][0] = F((-5)*i, [i-1, 0])
    return array


def compute(array, seq1, seq2, ms, mms, gs):
    
    """
    ms: Match Score
    mms: Mismatch Score
    gs: Gap Score
    
    """
    
    
    row, col = len(seq2), len(seq1)
    for i in range(1, row+1):
        for j in range(1, col+1):
            if seq1[j-1] == seq2[i-1]:
                s = ms
            else:
                s = mms  
            lu = [array[i-1][j-1].score + s, [i-1, j-1]]
            left = [array[i-1][j].score + gs, [i-1, j]]
            up = [array[i][j-1].score   + gs, [i, j-1]]
            max_choice = max([lu, left, up], key=lambda x: x[0])
            score= max_choice[0]
            pointer = max_choice[1]
            array[i][j] = F(score, pointer)
    return array


def backtrack(array, seq1, seq2):
    s1 = []
    s2 = []
    row, col = len(seq2), len(seq1)
#     while array[row][col].score != 0:
    while row > 0 or col > 0:
        i, j = array[row][col].pointer
        if i+1 == row and j+1 == col:
            s1.append(seq1[col-1])
            s2.append(seq2[row-1])
            row, col = i, j
        elif row == i+1 and col == j:
            s1.append("-")
            s2.append(seq2[i])
            row, col = i, j
        elif row == i and col == j+1:
            s1.append(seq1[j])
            s2.append("-")
            row, col = i, j
    s1 = ''.join(s1[::-1])
    s2 = ''.join(s2[::-1])
    return s1, s2


def align(seq1, seq2, ms=10, mms=-5, gs=-5):
    
    x, y = len(seq2)+1 , len(seq1)+1

    array = init_array(x, y)
    array = compute(array, seq1, seq2, ms, mms, gs)
    s1, s2 = backtrack(array, seq1, seq2)
    max_score = array[x-1][y-1].score
    
    return max_score, s1, s2
