from collections import Counter


def hamming(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Sequence lengths provided to hamming function are not equal.")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def qhamming(s1, s2, limit=1):
    if len(s1) != len(s2):
        raise ValueError("Sequence lengths provided to hamming function are not equal.")
    c = 0
    for ch1, ch2 in zip(s1, s2):
        if ch1 != ch2:
            c += 1
            if c > limit:
                return False
    return True


def determine_consensus(counts, threshold=0.7):
    cons = dict()
    for i, counter in enumerate(counts):
        c = counter.most_common(1)[0][1] / sum(counter.values())
        if c >= threshold:
            cons[i] = counter.most_common(1)[0][0]
    return cons


def gen_coord_map(s1, s2):
    coord_map = dict()
    ind_s2 = 0
    for i, s in enumerate(s1):
        if s == '-':
            coord_map[i] = None
            continue
        if s == s2[ind_s2]:
            coord_map[i] = ind_s2
            ind_s2 += 1
            continue
        if s != '-' and s2[ind_s2] == '-':
            while s2[ind_s2] == '-':
                ind_s2 += 1
            if s != s2[ind_s2]:
                raise Exception('Consensus sequences do not match!')
            else:
                coord_map[i] = ind_s2
                ind_s2 += 1
                continue
        raise Exception(f'Strange condition at ({i}, {ind_s2}), {s1[i-5:i+5]}, {s2[ind_s2-5:ind_s2+5]}')
    return coord_map


def get_counts(seq_dict):
    """
    Get counts of each nucleotide at each position (accepts gaps).

    Arguments:
        seq_dict (dict): Dictionary of sequences, should all be aligned and have equivalent lengths.
    Returns:
        List of Counter objects with length equal to the length of sequences.
    """
    ccs_len = len(next(iter(seq_dict.values())))
    return [Counter(''.join((seq_dict[seq][j] for seq in seq_dict))) for j in range(ccs_len)]
