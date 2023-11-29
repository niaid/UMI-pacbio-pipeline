import sys
import numpy as np

if __name__ == '__main__':
    stdin = sys.stdin
    counts = []
    for line in stdin:
        c = int(line.strip().split()[0])
        counts.append(c)
    counts = np.array(counts)
    print("\tSingleton UMIs: {} ({:1.1f}%)".format(len(np.argwhere(counts == 1)),
                                                   100 * len(np.argwhere(counts == 1)) / len(counts)))
    print("\tMean CCS reads per UMI:", '{:1.2f}'.format(np.mean(counts)))
    print("\tMedian CCS reads per UMI:", '{:1.1f}'.format(np.median(counts)))
    print("\tUMIs with 10 or more CCS:", '{} ({:1.1f}%)'.format(len(np.argwhere(counts >= 10)),
                                                                100 * len(np.argwhere(counts >= 10)) / len(counts)))
    print("\tLargest bin:", '{}'.format(np.max(counts)))
