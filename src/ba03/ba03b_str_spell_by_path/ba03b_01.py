'''
Rosalind: BA3B
String Spelled by a Genome Path Problem

Find the string spelled by a genome path.

Given: A sequence of k-mers Pattern1, ... , Patternn
such that the last k - 1 symbols of Patterni are equal
to the first k - 1 symbols of Patterni+1 for i from 1 to n-1.

Return: A string Text of length k+n-1 where the i-th k-mer in Text is equal to Patterni for all i.

Sample Dataset
ACCGA
CCGAA
CGAAG
GAAGC
AAGCT

Sample Output
ACCGAAGCT

-------------------------------------------------

solution by Edgar Martínez Encarnación
'''

import time

def main():
    f = open('/home/wsl/rosalind/data/ba03b.txt', 'r')
    reads=[]
    for line in f.readlines():
        reads.append(line.replace('\n',''))
    print(reads[0] + ''.join([item[-1] for item in reads[1:]]))
    f.close()

# main function
if __name__ == "__main__":
    main()