"""
Rosalind: BA3I
Find a k-Universal Circular String

A k-universal circular string is a circular string
that contains every possible k-mer constructed over a given alphabet.

k-Universal Circular String Problem
Find a k-universal circular binary string.

Given: An integer k.

Return: A k-universal circular string.
(If multiple answers exist, you may return any one.)

Sample Dataset
4

Sample Output
0000110010111101

═════════════════════════════════════════════════

https://github.com/wikiselev/coursera-bioinformatics-algorithms/blob/master/30-k-universal-circular-string-problem-cheated-with-python-script-from-wikipedia/30-k-universal-circular-string-problem.py
"""
#!/usr/bin/env python
import time


def de_bruijnn(k, n):
    """
    De Bruijn sequence for alphabet size k
    and subsequences of length n.
    """
    a = [0] * k * n
    sequence = []
    def db(t, p):
        if t > n:
            if n % p == 0:
                for j in range(1, p + 1):
                    sequence.append(a[j])
        else:
            a[t] = a[t - p]
            db(t + 1, p)
            for j in range(a[t - p] + 1, k):
                a[t] = j
                db(t + 1, t)
    db(1, 1)
    return sequence

seq = de_bruijnn(2, 17)
print(seq)


def main():
    f = open('/home/wsl/rosalind/data/ba03i.txt', 'r')
    k = int(f.readline().strip())
    f.close()

    start_time = time.time()
    print(k_universal_strings(2))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
