'''
Compare two lists with duplicates

Python | Difference of two lists including duplicates
https://www.geeksforgeeks.org/python-difference-of-two-lists-including-duplicates/
'''

# 1.

from collections import Counter

test_list1 = [1, 3, 4, 5, 1, 3, 3]
test_list2 = [1, 3, 5]
print("The original list 1 : " + str(test_list1))
print("The original list 2 : " + str(test_list2))

res = list((Counter(test_list1) - Counter(test_list2)).elements())
print("The list after performing the subtraction : " + str(res))

# 2.

test_list1 = [1, 3, 4, 5, 1, 3, 3]
test_list2 = [1, 3, 5]
print("The original list 1 : " + str(test_list1))
print("The original list 2 : " + str(test_list2))

res = map(lambda x: test_list1.remove(x) if x in test_list1 else None, test_list2)
print("The list after performing the subtraction : " + str(test_list1))
