"""
Rosalind: BA5C
Find a Longest Common Subsequence of Two Strings

Longest Common Subsequence Problem
Given: Two strings.

Return: A longest common subsequence of these strings.

Sample Dataset
AACCTTGG
ACACTGTGA

Sample Output
AACTGG

═════════════════════════════════════════════════

solution by Matthew Burch
"""

#!/usr/bin/env python3

import time

def LongestCommonSubsequence(str1, str2):
    """
    Inspired by:  Tushar Roy's Dynamic Programming LCS video
    https://www.youtube.com/watch?v=NnD96abizww
    """
    l1, l2 = len(str1), len(str2)
    arr = [i[:] for i in [[0] * (l2+1)] * (l1+1)]  # http://stackoverflow.com/a/2398187/4427457
    for idx1,s1 in enumerate(str1):
        for idx2,s2 in enumerate(str2):
            if s1 == s2:
                arr[idx1+1][idx2+1] = arr[idx1][idx2] + 1
            else:
                arr[idx1+1][idx2+1] = max(arr[idx1][idx2+1], arr[idx1+1][idx2])
    subsequence = []
    i,j = l1,l2
    while i > 0 and j > 0:
        curr = arr[i][j]
        up = arr[i-1][j]
        left = arr[i][j-1]
        if curr > left and curr > up:
            subsequence.insert(0, str1[i-1])
            i,j = i-1,j-1
        elif curr > left:
            i -= 1
        else:
            j -= 1
    return "".join(subsequence)


# main function
def main():
    try:
        with open('/home/wsl/rosalind/data/ba05c.txt', 'r') as f:
            lines = [line.strip() for line in f.readlines()]
            str1 = lines[0]
            str2 = lines[1]
    except OSError as e:
        print("!!!!! No such file or directory. !!!!!")
        print(e.errno)

    start_time = time.time()
    print(LongestCommonSubsequence(str1, str2))
    print("--- %s seconds ---" % (time.time() - start_time))


# execute main fx
if __name__ == "__main__":
    main()