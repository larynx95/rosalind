"""
Rosalind: BA5A
Find the Minimum Number of Coins Needed to Make Change

The Change Problem
Find the minimum number of coins needed to make change.

Given: An integer money and an array Coins of positive integers.

Return: The minimum number of coins with denominations Coins that changes money.

Sample Dataset
40
1,5,10,20,25,50

Sample Output
2

═════════════════════════════════════════════════

    [ Where am I? ]

    * HERE: Dynamic programming (DP) in "coin change problem" (BA5A)
      ↓
    * NEXT: Applying dynamic programming to "Manhattan tourist problem" - Longest Path (BA5B)

Plan 1.
- Can I use one coin twice or more?  <-- This!
  Can I use one coin only once?      <-- Not this!
- GreedyChange algorithm
  Pseudoce (textbook, incorrect):
  ╔══════════════════════════════════════════════════════════════════════════════╗
  ║ GREEDYCHANGE(money)                                                          ║
  ║     Change <- empty collection of coins                                      ║
  ║     while money > 0                                                          ║
  ║         coin <- largest denomination that is less than or equal to money     ║
  ║         add a coin with denomination coin to the collection of coins Change  ║
  ║         money <- money - coin                                                ║
  ║     return Change                                                            ║
  ╚══════════════════════════════════════════════════════════════════════════════╝

- dividing money by the largest coin value, again and again

  40    1,5,10,20,25,50
  40 -> [25] * 1 + 15
  15 -> [10] * 1 + 5
  5  -> [ 5] * 1 + 0       <-- what if not zero?

- find the largest coin value that is less than or equal to money
  i. simply, the last element of sorted list of coins
  ii. binary search (not good, TODO: try this later)

- Is this 'greedyChange' approach always correct? No.

  48 = 40 + 5 + 1 + 1 + 1  <-- incorect
     = 24 + 24             <-- correct

Plan 2.
- another GreedyChange approach, but different
  get ALL possible cases, then find a case with minimum count
- simulation
  48, [120,40,30,24,20,10, 5, 4, 1]    columns = len(coins)
        0   0  0  0  0  0  0  0 48     rows    = ??
        0   0  0  0  0  0  0  1 44
        0   0  0  0  0  0  0  2 40     i. for-loops
                              ...      ii. recursion
        0   0  0  0  0  0  0 12  0     iii. some other way
                           ...
        0   0  0  0  0  0  9 ...
                        ...
        0   0  0  0  0  4 ...
                     ...
        0   0  0  0  2 ...
                  ...
        0   0  0  2 ...
        0   0  1 ...
        0   1 ...
- TODO: How can I get all the possible cases?

Plan 3.
- recursive algorithm

  part A                                             part B (important)
  -------------------------------------------------------------------------
  a minimal collection of coins totaling 75 denarii, plus a 1-denarius coin
  a minimal collection of coins totaling 72 denarii, plus a 4-denarius coin
  a minimal collection of coins totaling 71 denarii, plus a 5-denarius coin

                           ┌ MinNumCoins(money - coin_1) + 1
  MinNumCoins(money) = min ┤
                           └ MinNumCoins(money - coin_d) + 1

- Pseudocode (textbook):
  This function is not useful when money is large value.
  ╔══════════════════════════════════════════════════════════════════╗
  ║ RECURSIVECHANGE(money, COINS)                                    ║
  ║     if money = 0                                                 ║
  ║         return 0                                                 ║
  ║     minNumCoins <- infinite                                      ║
  ║     for i <- 1 to |COINS|                                        ║
  ║         if money >= coin_i                                       ║
  ║             numCoins <- RECURSIVECHANGE(money - coin_i, COINS)   ║
  ║             if numCoins + 1 < minNumCoins                        ║
  ║                 minNumCoins <- numCoins + 1                      ║
  ║     return minNumCoins                                           ║
  ╚══════════════════════════════════════════════════════════════════╝

Plan 4.
- dynamic programming (calculate once)
  coins = [1,4,5], money = [1..12]

                           ┌ MinNumCoins(m - 5) + 1
  MinNumCoins(money) = min ┼ MinNumCoins(m - 4) + 1
                           └ MinNumCoins(m - 1) + 1

  MinNumCoins(0)     =                             => 0      []

  MinNumCoins(1)     =       MinNumCoins(1-1) + 1  => 1      [1]

  MinNumCoins(2)     =       MinNumCoins(2-1) + 1  => 2      [1,1]

  MinNumCoins(3)     =       MinNumCoins(3-1) + 1  => 3      [1,1,1]

  MinNumCoins(4)     = min ┌ MinNumCoins(4-4) + 1  => 1      [4]
                           └ MinNumCoins(4-1) + 1  => 4      [1,1,1,1]

                           ┌ MinNumCoins(5-5) + 1  => 1      [5]
  MinNumCoins(5)     = min ┼ MinNumCoins(5-4) + 1  => 2      [1,4]
                           └ MinNumCoins(5-1) + 1  => 5      [1,1,1,1,1]

                           ┌ MinNumCoins(6-5) + 1  => 2      [1,5]
  MinNumCoins(6)     = min ┼ MinNumCoins(6-4) + 1  => 2      [1,1,4]
                           └ MinNumCoins(6-1) + 1  => 2      [5,1]

                           ┌ MinNumCoins(7-5) + 1  => 3      [1,1,5]
  MinNumCoins(7)     = min ┼ MinNumCoins(7-4) + 1  => 4      [1,1,1,4]
                           └ MinNumCoins(7-1) + 1  => 3      [5,1,1]

              m : 0   1   2   3   4   5   6   7   8   9   10   11   12
  MinNumCoins(m): 0   1   2   3   1   1   2   3   2   2    2    3    3

- Pseudocode (textbook):
  ╔════════════════════════════════════════════════════════════════════╗
  ║ DPCHANGE(money, COINS)                                             ║
  ║     MINNUMCOINS(0) <- 0                                            ║
  ║     for m <- 1 to money                                            ║
  ║         MINNUMCOINS(m) 1                                           ║
  ║         for i <- 1 to |COINS|                                      ║
  ║             if m >= coin_i                                         ║
  ║                 if MINNUMCOINS(m - coin_i) + 1 < MINNUMCOINS(m)    ║
  ║                     MINNUMCOINS(m) <- MINNUMCOINS(m - coin_i) + 1  ║
  ║     return MINNUMCOINS(money)                                      ║
  ╚════════════════════════════════════════════════════════════════════╝

- TODO: (question in textbook)
  STOP and Think: If money = 10^9,
  DPCHANGE requires a huge array of size 10^9.
  Modify the DPCHANGE algorithm
  so that the array size required does not exceed the value of the largest coin denomination.

- [DONE]: (question in textbook)
  STOP and Think: Recall that our original goal was to make change,
  not just compute MINNUMCOINS(money).
  Modify DPCHANGE so that it not only computes the minimum number of coins
  but also returns these coins.

═════════════════════════════════════════════════

References:
- How to find the highest number less than target value in a list?
  https://stackoverflow.com/questions/20023004/how-to-find-the-highest-number-less-than-target-value-in-a-list
- Binary Search (bisect) in Python
  https://www.geeksforgeeks.org/binary-search-bisect-in-python/
- Find the smallest number that is greater than a given number in a sorted list
  https://stackoverflow.com/questions/13669770/find-the-smallest-number-that-is-greater-than-a-given-number-in-a-sorted-list
- Recursive generator for change money - python
  https://stackoverflow.com/questions/44575298/recursive-generator-for-change-money-python
- Gist @DiegoGallegos4/change.py
  https://gist.github.com/DiegoGallegos4/9f9e1090f95ec27963fc8088b483e9a6
"""

import time


#################################################
# change function
#################################################


def change_greedy(money, rsorted_coins):
    """
    (int,[int]) -> [int]
    algorithm in textbook (incorrect result)
    rsorted_coins: a list of reversely sorted integers
    >>> change_greedy(40,[50,25,20,10,5,1])
        [25,10,5]
    >>> change_greedy(48,[120,40,30,24,20,10,5,4,1])
        [40,5,1,1,1]    # <-- This is not the fewest number of coins. [24,24]
    """
    change = []
    while money > 0:
        # get the largest coin that is less than or equal to money
        largest_coin = 0
        for coin in rsorted_coins:
            if coin <= money:    # <-- '<' makes infinite loop, should be '<='
                largest_coin = coin
                break
        # update the value of money and append it to the result list
        money -= largest_coin
        change.append(largest_coin)
    return change


def change_recursive(money, coins):
    """
    (int,[int]) -> [int]
    returns a list of changes
    >>> change_recursive(11,[5,3,1])
        [5,5,1]
    """
    if money == 0:
        return []
    min_changes = [min(coins)] * money
    for i in range(len(coins)):
        if money >= coins[i]:
            changes = change_recursive(money - coins[i], coins)
            if len(changes) + 1 < len(min_changes):
                min_changes = [coins[i]] + changes
    return min_changes


def change_dynamic(money, coins):
    """
    (int,[int]) -> [int]
    >>> change_dynamic(11,[5,3,1])
        [1, 5, 5]
    >>> change_dynamic(48,[120,40,30,24,20,10,5,4,1])
        [24, 24]
    """
    min_num_coins = [0] * (money + 1)    # index as money, elem as MinNumCoins
    min_ls_coins = [[]] * (money + 1)    # index as money, elem as list of coins with MinNumCoins
    for m in range(1, money+1):
        min_num_coins[m] = float('inf')
        for i in range(len(coins)):
            if m >= coins[i]:
                if min_num_coins[m - coins[i]] + 1 < min_num_coins[m]:
                    min_num_coins[m] = min_num_coins[m - coins[i]] + 1  # swap minimum number of coins
                    min_ls_coins[m]  = min_ls_coins[m - coins[i]] + [coins[i]] # swap list of coins with minimum number of coins
    return min_ls_coins[money]


def change_dynamic_better(money, coins):
    """
    (int,[int]) -> [int]
    TODO: I haven't fully understood dynamic programming yet.
    >>> change_dynamic_better(11,[5,3,1])
        [1, 5, 5]
    >>> change_dynamic_better(48,[120,40,30,24,20,10,5,4,1])
        [24, 24]
    """
    min_ls_coins = [[]] * (money + 1)  # <-- TODO: too many! Fix this!
    for m in range(1, money+1):
        min_num_coins = float('inf')
        for i in range(len(coins)):
            if m >= coins[i]:
                temp_ls = min_ls_coins[m - coins[i]] + [coins[i]]
                if len(temp_ls) < min_num_coins:  # <-- use 'len', not 'sum'
                    min_ls_coins[m]  = temp_ls
                    min_num_coins = len(temp_ls)
    return min_ls_coins[money]


#################################################
# minimum number of changes
#################################################


def change_greedy_recursive(money, coins):
    """
    (int,[int]) -> [int]
    recursive version of GreedyChange
    returns incorrect result too, TODO: Fix this.
    >>> change_greedy_recursive(48,[120,40,30,24,20,10,5,4,1])
        [0,1,0,0,0,0,1,0,3]    # <-- incorect, it should be [0,0,0,2,0,0,0,0,0]
    """
    if money <= 0:         # base case, if no money
        return []
    if len(coins) == 0:    # base case, if no coin
        return []
    num = money // coins[0]
    if num == 0:           # recursive case, if the largest coin > money
        return [0] + change_greedy_recursive(money, coins[1:])
    else:                  # recursive case, if the largest coin <= money
        return [num] + change_greedy_recursive(money-coins[0]*num, coins[1:])


def min_num_change_recursive(money, coins):
    """
    (int,[int]) -> int
    algorithm is textbook
    >>> min_num_change_recursive(11,[5,3,1])
        3
    >>> min_num_change_recursive(40,[50,25,20,10,5,1])
        2
    >>> min_num_change_recursive(48,[120,40,30,24,20,10,5,4,1])
        TODO: Fix this. It takes too long time! Not working. Why?
    """
    if money == 0:
        return 0
    min_ncoins = float('inf')
    for i in range(len(coins)):
        if money >= coins[i]:
            ncoins = min_num_change_recursive(money-coins[i], coins)
            if ncoins + 1 < min_ncoins:
                min_ncoins = ncoins + 1
    return min_ncoins


def min_num_change_dynamic(money, coins):
    """
    (int,[int]) -> int
    dynamic programming algorithm in textbook
    TODO: This is awesome! Find more resources about dynamic programming.
    >>> min_num_change_dynamic(48,[120,40,30,24,20,10,5,4,1])
        2
    """
    min_num_coins = [0] * (money + 1)
    for m in range(1, money+1):
        min_num_coins[m] = float('inf')
        for coin in coins:
            if m >= coin:
                if min_num_coins[m - coin] + 1 < min_num_coins[m]:
                    min_num_coins[m] = min_num_coins[m - coin] + 1
    return min_num_coins[money]


def main():
    try:
        with open('/home/wsl/rosalind/data/ba05a.txt', 'r') as f:
            lines = [line.strip() for line in f.readlines()]
            money = int(lines[0])
            coins = [int(line.strip()) for line in lines[1].split(',')]
    except OSError as e:
        print("!!!!! No such file or directory. !!!!!")
        print(e.errno)

    start_time = time.time()
    print(min_num_change_dynamic(money, coins))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
