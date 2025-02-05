#  Copyright (c) National Institutes of Health
#  Written by Pierce Radecki

import itertools
import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import networkx
from Bio import pairwise2
import seqlib
import Levenshtein
from Bio.SeqIO.FastaIO import SimpleFastaParser
from itertools import combinations
from collections import Counter
import scipy.stats


def determine_consensus(counts, threshold=0.7):
    cons = dict()
    for i, counter in enumerate(counts):
        c = counter.most_common(1)[0][1] / sum(counter.values())
        if c >= threshold:
            cons[i] = counter.most_common(1)[0][0]
    return cons


def get_dm_sites(counts):
    dms = []
    for i, c in enumerate(counts):
        for letter in c:
            if c[letter] > 1 and c[letter] <= min((sum(c.values())*0.5, max((0.01*sum(c.values()), 10)))):
                dms.append((i, letter))
    return dms


if __name__ == "__main__":

    project = sys.argv[1]
    sample = sys.argv[2]

    with open(os.path.join(project, f'sgs-prelim/mafft/{sample}.fasta'), 'r') as f:
        sgs_seqs = {record[0]: record[1] for record in SimpleFastaParser(f)}

    sgs_len = len(next(iter(sgs_seqs.values()))) - 10
    sgs_counts = [Counter(''.join((sgs_seqs[s][j] for s in sgs_seqs))) for j in range(sgs_len)]
    sgs_cons = determine_consensus(sgs_counts, threshold=0.0)

    dm_sites = get_dm_sites(sgs_counts)

    # sys.exit()
    data_umi = dict()

    for sgs_umi in sgs_seqs:

        sgs_name_umi = sgs_umi
        umi = sgs_umi.split('_')[-2]
        umi_size = int(sgs_umi.split('_')[-1])

        muts = []
        shift = 0
        for i in sgs_cons:
            si = sgs_seqs[sgs_name_umi][i]
            sc = sgs_cons[i]
            if sc == '-':
                shift += 1
            if si != sc:
                muts.append((i - shift, si, sc))

        muts_consolidated_prelim = []

        open_del = False
        del_start = None
        del_seq = ''

        open_ins = False
        ins_start = None
        ins_seq = ''

        for mut in muts:
            if mut[1] == '-' and mut[2] != '-':
                if not open_del:
                    open_del = True
                    del_start = mut[0]
                    del_seq = mut[2]
                else:
                    del_seq += mut[2].replace('-', '')
            elif mut[1] != '-' and mut[2] == '-':
                if not open_ins:
                    open_ins = True
                    ins_start = mut[0]
                    ins_seq = mut[1]
                else:
                    ins_seq += mut[1]
            else:
                if open_del:
                    muts_consolidated_prelim.append(('del', del_start, del_seq, len(del_seq)))
                    open_del = False
                elif open_ins:
                    muts_consolidated_prelim.append(('ins', ins_start, ins_seq, len(ins_seq)))
                    open_ins = False

            if mut[1] not in '-' and mut[2] not in '-':
                muts_consolidated_prelim.append(('snv', mut[0], mut[1], mut[2]))

        if open_del:
            muts_consolidated_prelim.append(('del', del_start, del_seq, len(del_seq)))
            open_del = False
        elif open_ins:
            muts_consolidated_prelim.append(('ins', ins_start, ins_seq, len(ins_seq)))
            open_ins = False

        muts_consolidated = []
        muts_consolidated_not_dm = []
        for mut in muts_consolidated_prelim:
            if mut[0] == 'snv':
                if (mut[1], mut[2]) in dm_sites:
                    muts_consolidated.append(mut)
                else:
                    muts_consolidated_not_dm.append(mut)
                    continue
            if mut[0] == 'del':
                if (mut[1], '-') in dm_sites:
                    muts_consolidated.append(mut)
                else:
                    muts_consolidated_not_dm.append(mut)
                    continue

        fp_parent_align = os.path.join(project, f'bins/{sample}/aligned/{umi}-aligned.fasta')
        mfs = []

        with open(fp_parent_align, 'r') as f:
            sfp = SimpleFastaParser(f)
            cons = next(sfp)[1]
            seqs = {record[0]: record[1] for record in sfp}

        map = pairwise2.align.globalms(cons, sgs_seqs[sgs_name_umi].replace('-', ''), 2, -1, -5, -.6, one_alignment_only=True,
                                       penalize_end_gaps=False)[0]

        map_sgs2bin = dict()
        for i, (a, b) in enumerate(zip(map.seqB, map.seqA)):
            if a == b:
                if b == '-':
                    continue
                map_sgs2bin[i] = len(map.seqB[:i].replace('-', ''))
                continue

        sgs_len = len(next(iter(seqs.values())))
        counts = [Counter(''.join((seqs[s][j] for s in seqs))) for j in range(sgs_len)]

        data_errs = []

        for i in map_sgs2bin:
            c = counts[map_sgs2bin[i]]
            mc = c.most_common()[0][0].lower()

            errs = [c[base] for base in ('a', 'c', 'g', 't') if base != mc and base in c]
            if errs:
                data_errs.append(max(errs))
            else:
                continue

        if data_errs:
            thresh = max([1.0, np.median(data_errs), 0.001 * len(seqs)]) + 1.0
        else:
            thresh = 2.0

        for i in map_sgs2bin:
            c = counts[map_sgs2bin[i]]
            if len(c) == 1 or (len(c) == 2 and '-' in c):
                continue
            else:
                for base in filter(lambda a: a not in (cons[map_sgs2bin[i]], '-'), c.keys()):
                    mf = c[base]
                    if mf > thresh:
                        mfs.append((i, base))

        data_umi[umi] = [umi_size, muts_consolidated, mfs, muts_consolidated_not_dm]

    G = networkx.Graph()
    G.add_nodes_from(data_umi)

    Gfp = networkx.Graph()
    Gfp.add_nodes_from(data_umi)

    fb_ed1 = set()
    fp_dict = {umi: set() for umi in data_umi}

    for umi1, umi2 in itertools.combinations(data_umi, r=2):
        if umi1 == umi2:
            continue
        if Levenshtein.distance(umi1, umi2) == 1:
            if data_umi[umi1][0] > 2 * data_umi[umi2][0] or data_umi[umi1][0] < 0.5 * data_umi[umi2][0]:
                G.add_edge(umi1, umi2)
                smaller_bin = min([umi1, umi2], key=lambda a: data_umi[a][0])
                larger_bin = max([umi1, umi2], key=lambda a: data_umi[a][0])
                fb_ed1.add(smaller_bin)
                fp_dict[smaller_bin].add(larger_bin)
                print(smaller_bin, data_umi[smaller_bin][0], '-->', larger_bin, data_umi[larger_bin][0])


    def match(m1, m2):
        # Method to determine if two bins are likely associated as UMI bins from the same RNA
        # The larger bin will be kept and called as the true bin
        larger_bin = max([m1, m2], key=lambda a: a[0])
        smaller_bin = min([m1, m2], key=lambda a: a[0])
        bin1_rt = [(m[1], m[2]) for m in larger_bin[1]]
        bin1_pcr = larger_bin[2]
        bin1_aux = [(m[1], m[2]) for m in larger_bin[3]]
        bin2_rt = [(m[1], m[2]) for m in smaller_bin[1]]
        bin2_aux = [(m[1], m[2]) for m in smaller_bin[3]]
        if len(bin2_rt) >= 1 or len(bin2_aux) >= 1:
            if all([(m in bin1_rt) or (m in bin1_pcr) for m in bin2_rt]) and all([m in bin2_rt for m in bin1_rt]):
                if all([(m in bin1_aux) or (m in bin1_pcr) for m in bin2_aux]):
                    return True
        return False


    for umi1, umi2 in itertools.combinations(data_umi, r=2):

        if umi1 == umi2:
            continue
        sgs_name_umi1 = next(filter(lambda a: a.__contains__(f'_{umi1}'), sgs_seqs.keys()))
        sgs_name_umi2 = next(filter(lambda a: a.__contains__(f'_{umi2}'), sgs_seqs.keys()))

        if Levenshtein.distance(umi1, umi2) >= 3:
            if (umi1 in fb_ed1) or (umi2 in fb_ed1):
                continue

        if sum([a != b and '-' not in (a, b) for a, b in zip(sgs_seqs[sgs_name_umi1], sgs_seqs[sgs_name_umi2])]) <= 5 and (
                data_umi[umi1][0] > 10 * data_umi[umi2][0] or data_umi[umi1][0] < 0.1 * data_umi[umi2][0]):
            if match(data_umi[umi1], data_umi[umi2]):
                smaller_bin = min([umi1, umi2], key=lambda a: data_umi[a][0])
                larger_bin = max([umi1, umi2], key=lambda a: data_umi[a][0])
                fb_ed1.add(smaller_bin)
                fp_dict[smaller_bin].add(larger_bin)
                Gfp.add_edge(umi1, umi2)
                print(smaller_bin, data_umi[smaller_bin], '-->', larger_bin, data_umi[larger_bin])

    fake_bins = set()
    for umi in fp_dict:
        if len(fp_dict[umi]) > 1:
            fake_bins.add(umi)
            G.remove_node(umi)
            Gfp.remove_node(umi)

    # Generate figure showing association via edit distance and PCR fingerprinting
    f = plt.figure()

    ccs = list(networkx.connected_components(G))

    grid_size = np.round(np.sqrt(len(ccs)))

    pos_all = dict()
    for i, cc in enumerate(ccs):
        pos = networkx.spring_layout(cc, center=(i // grid_size, i % grid_size), scale=0.1, iterations=1000)
        pos_all.update(pos)
        Gcc = networkx.subgraph_view(G, filter_node=lambda a: a in cc)
        networkx.draw_networkx(Gcc, pos=pos, node_size=[10 * np.sqrt(data_umi[d][0]) for d in Gcc.nodes], font_size=5,
                               verticalalignment='top',
                               node_color=['blue' if data_umi[d][0] > 1000 else 'red' for d in Gcc.nodes],
                               font_color='black')

    networkx.draw_networkx_edges(Gfp, pos=pos_all, edge_color='magenta')
    plt.savefig(os.path.join(project, f'sgs-prelim/fake-umi/{sample}-fp.png'), format='png', dpi=1000)
    plt.close(f)

    G_new = networkx.Graph()
    G_new.add_nodes_from(G)
    G_new.add_edges_from([(a[0], a[1]) for a in G.edges])
    G_new.add_edges_from([(a[0], a[1]) for a in Gfp.edges])

    ccs = list(networkx.connected_components(G_new))
    grid_size = np.round(np.sqrt(len(ccs)))

    pos_all = dict()
    real_bins = [max(data_umi, key=lambda umi: data_umi[umi][0])]
    fake_bins = [umi for umi in fake_bins]

    for i, cc in enumerate(ccs):
        pos = networkx.spring_layout(cc, center=(i // grid_size, i % grid_size), scale=0.1, iterations=1000)
        pos_all.update(pos)
        real = max(cc, key=lambda a: data_umi[a][0])
        real_bins.append(real)
        for cci in cc:
            if cci == real:
                continue
            else:
                fake_bins.append(cci)

    f = plt.figure()

    # Perform filtering based on mixture model of fake/real bins
    # Initialize model with known values
    real_bins_iter = [umi for umi in data_umi if umi not in fake_bins]

    real_bin_sizes = [data_umi[umi][0] for umi in data_umi if umi in real_bins]
    fake_bin_sizes = [data_umi[umi][0] for umi in data_umi if umi in fake_bins] + [9]
    # Add 9 as placeholder representing smaller fake bins

    fake_bins.append('X')
    U = [umi for umi in data_umi] + ['X']
    X = np.array([data_umi[umi][0] for umi in data_umi] + [9])

    mu = np.mean(real_bin_sizes)
    dr = scipy.stats.norm(loc=mu, scale=np.std(real_bin_sizes))
    df = scipy.stats.expon(scale=np.mean(fake_bin_sizes))
    x = np.arange(1, max(real_bin_sizes))
    plt.hist(X, color='grey', alpha=0.7)

    P = np.array([dr.pdf(x) / (df.pdf(x) + dr.pdf(x)) for x in X])
    for i, umi in enumerate(U):
        if umi in fake_bins:
            P[i] = 0
        else:
            continue

    hist_bins_size = (max(X) - min(X)) / 10
    plt.plot(x, dr.pdf(x) * hist_bins_size * sum(p for p in P), color='g', ls='--', lw=2)
    plt.plot(x, df.pdf(x) * hist_bins_size * sum(1-p for p in P), color='maroon', ls='--', lw=2)

    tw2 = plt.gca().twinx()
    tw2.plot(x, dr.pdf(x)*sum(p for p in P) / (df.pdf(x)*sum(1-p for p in P) + dr.pdf(x)*sum(p for p in P)),
             color='k')
    tw2.set_ylim((-0.05, 1.05))
    tw2.set_ylabel('Prob(real)')
    # Initial model snapshot
    plt.savefig(os.path.join(project, f'sgs-prelim/fake-umi/{sample}-init-em.pdf'), format='pdf')

    lls = [1]

    ll = sum([np.log(p*dr.pdf(x) + (1-p)*df.pdf(x)) for p, x in zip(P, X)])
    lls.append(ll)

    iteration_count = 1
    while np.abs((lls[-1] - lls[-2]) / lls[-2]) > 0.01:

        if iteration_count > 100:
            # Accept converged solution after 100 iterations
            print('Maximum number of iterations reached.')
            break

        f = plt.figure()

        mean_pseudo_real = sum((p*x) for p, x in zip(P, X)) / sum(p for p in P)
        sigma_pseudo_real = np.sqrt(sum(p*(x-mu)**2 for p, x in zip(P, X)) / sum(p for p in P))

        # Constrain estimation of std. dev. / variance
        if sigma_pseudo_real < 1:
            sigma_pseudo_real = 1

        mu = mean_pseudo_real

        mean_pseudo_fake = sum((1-p)*x for p, x in zip(P, X)) / sum(1-p for p in P)

        dr = scipy.stats.norm(loc=mean_pseudo_real, scale=sigma_pseudo_real)
        df = scipy.stats.expon(scale=mean_pseudo_fake)
        x = np.arange(1, max(X))
        plt.hist(X, color='grey', alpha=0.7)

        hist = (max(X) - min(X)) / 10
        plt.plot(x, dr.pdf(x) * hist * sum(p for p in P), color='g', ls='--', lw=2)
        plt.plot(x, df.pdf(x) * hist * sum(1-p for p in P), color='maroon', ls='--', lw=2)

        tw2 = plt.gca().twinx()
        tw2.plot(x, dr.pdf(x)*sum(p for p in P) / (df.pdf(x)*sum(1-p for p in P) + dr.pdf(x)*sum(p for p in P)), color='k')
        tw2.set_ylim((-0.05, 1.05))
        tw2.set_ylabel('Prob(real)')

        P = np.array([dr.pdf(x)*sum(p for p in P) / (df.pdf(x)*sum(1-p for p in P) + dr.pdf(x)*sum(p for p in P)) for x in X])
        for i, umi in enumerate(U):
            # Enforce hard constraints of UMIs known to be fake
            if X[i] > dr.mean():
                P[i] = 1
            if umi in fake_bins:
                P[i] = 0
            else:
                continue

        ll = sum([np.log(p*dr.pdf(x) + (1-p)*df.pdf(x)) for p, x in zip(P, X)])
        lls.append(ll)

        real_bins_iter = []
        for i, x in enumerate(X):
            if U[i] == 'X':
                continue
            if P[i] >= 0.5:  # Threshold of 0.5
                print('Real bin:', U[i], data_umi[U[i]][0], P[i])
                real_bins_iter.append(U[i])
            else:
                print('Fake bin:', U[i], data_umi[U[i]][0], P[i])

        # Save snapshot of model at iteration
        plt.savefig(os.path.join(project, f'sgs-prelim/fake-umi/{sample}-{iteration_count}-em.pdf'), format='pdf')
        iteration_count += 1
        plt.close(f)

    # Save final sequences that are predicted likely real
    with open(os.path.join(project, f'sgs/{sample}_final.fasta'), 'w') as f:
        for umi in sorted(real_bins_iter, key=lambda a: data_umi[a][0], reverse=True):
            for seq in sgs_seqs:
                if umi in seq:
                    s = sgs_seqs[seq].replace('-', '')
                    f.write(f'>{seq}\n{s}\n')
