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

Plan 1.
- Can I use one coin twice or more?  <-- This!
  Can I use one coin only once?      <-- Not this!
- GreedyChange algorithm

    Pseudoce (textbook):
    GREEDYCHANGE(money)
        Change <- empty collection of coins
        while money > 0
            coin <- largest denomination that is less than or equal to money
            add a coin with denomination coin to the collection of coins Change
            money <- money - coin
        return Change

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
- algorithms in textbook
  ╔══════════════════════════════════════════════════════════════════════════════╗
  ║ GREEDYCHANGE(money)                                                          ║
  ║     Change <- empty collection of coins                                      ║
  ║     while money > 0                                                          ║
  ║         coin <- largest denomination that is less than or equal to money     ║
  ║         add a coin with denomination coin to the collection of coins Change  ║
  ║         money <- money - coin                                                ║
  ║     return Change                                                            ║
  ╚══════════════════════════════════════════════════════════════════════════════╝

  ╔══════════════════════════════════════════════════════════════════╗
  ║  RECURSIVECHANGE(money, COINS)                                   ║
  ║      if money = 0                                                ║
  ║          return 0                                                ║
  ║      minNumCoins <- infinite                                     ║
  ║      for i <- 1 to |COINS|                                       ║
  ║          if money >= coin_i                                      ║
  ║              numCoins <- RECURSIVECHANGE(money - coin_i, COINS)  ║
  ║              if numCoins + 1 < minNumCoins                       ║
  ║                  minNumCoins <- numCoins + 1                     ║
  ║      return minNumCoins                                          ║
  ╚══════════════════════════════════════════════════════════════════╝

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

═════════════════════════════════════════════════

- solution by Gist @DiegoGallegos4/change.py
  https://gist.github.com/DiegoGallegos4/9f9e1090f95ec27963fc8088b483e9a6
"""

import time


def find_max_coin(money, coins):
    max_coin = float('-inf')
    for coin in coins:
        if money - coin >= 0 and coin > max_coin:
            max_coin = coin
    return max_coin


def change_greedy(money, coins):
    """
    Running time:
    O(mn) where n is the number of coins and m is range of 0 to money
    """
    change = []
    while money > 0:
        coin = find_max_coin(money, coins)
        change.append(coin)
        money = money - coin
    return change


def change_greedy_sorted(money, coins):
    """
    Running time:
    O(max(n, m))
    """
    change = []
    sorted_coins = sorted(coins, reverse=True)
    for coin in coins:
        while money - coin >= 0:
            money -= coin
            change.append(coin)
    return coin


def change_recursive(money, coins):
    """
    Find minimum number of coins to give change for money
    Args:
    (array) coins denominations
    (float) money for which change is needed

    Returns:
    Minimum number of coins

    Running time:
    O(m*n^2)
    """
    if money == 0:
        return 0
    min_num_coins = float('inf')
    for i in range(len(coins)):
        if money >= coins[i]:
            num_coins = change_recursive(money - coins[i], coins)
            if num_coins + 1 <= min_num_coins:
                min_num_coins = num_coins + 1
    return min_num_coins


def change_dp(money, coins):
    """
    Dynamic Programming
    Running time:
    O(mn)
    """
    min_num_coins = [0] * (money + 1)
    for m in range(1, money + 1):
        min_num_coins[m] = float('inf')
        for coin in coins:
            if m >= coin:
                num_coins = min_num_coins[m - coin] + 1
                if num_coins < min_num_coins[m]:
                    min_num_coins[m] = num_coins
    return min_num_coins, min_num_coins[money]


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
    print(change_dp(money, coins))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
