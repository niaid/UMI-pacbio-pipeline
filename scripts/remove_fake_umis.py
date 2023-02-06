import json
import sys
import os
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
import seqlib
from Bio.SeqIO.FastaIO import SimpleFastaParser
from itertools import combinations


def check_fingerprint_match(umi1, umi2, fpdict, glmap, lgmap):
    fp1 = fpdict[umi1]['fp_pcr']
    fp2 = fpdict[umi2]['fp_pcr']
    dd1 = fpdict[umi1]['fp_dd']
    dd2 = fpdict[umi2]['fp_dd']

    fp2_local_ind = []
    fp2_global_ind = []
    for fp in fp2:
        fp2_local_ind_i = []
        fp2_global_ind_i = []
        for fpi in fp:
            g_ind = lgmap[umi2][str(fpi[0])]
            if g_ind:
                fp2_global_ind_i.append((g_ind, fpi[1]))
                fp1_ind = glmap[umi1][str(g_ind)]
                if fp1_ind:
                    fp2_local_ind_i.append((fp1_ind, fpi[1]))
            else:
                break
        if fp2_local_ind_i:
            fp2_local_ind.append(fp2_local_ind_i)
        if fp2_global_ind_i:
            fp2_global_ind.append(fp2_global_ind_i)

    for fp in fp1:
        for fp1i in fp:
            for fp2i in fp2_local_ind:
                if fp1i[0] == fp2i[0]:
                    print(umi1, umi2, 'PCR-PCRf', fp1i, '-->', fp2i)
                    return True

    for dd in dd1:
        for fp2i in fp2_global_ind:
            for fp2ij in fp2i:
                if dd[0] == fp2ij[0]:
                    print(umi1, umi2, 'dd-PCRf', dd1, '-->', fp2i)
                    return True

    for dd in dd1:
        for dd2i in dd2:
            if dd[0] == dd2i[0]:
                print(umi1, umi2, 'dd-PCRf', dd1, '-->', dd2)
                return True


if __name__ == "__main__":

    project = sys.argv[1]
    sample = sys.argv[2]

    with open(os.path.join(project, f'sgs-prelim/mafft/{sample}.fasta'), 'r') as f:
        sgs_seqs = {'_'.join(record[0].split('_')[-2:]): str(record[1]) for record in SimpleFastaParser(f)}

    with open(os.path.join(project, f'sgs-prelim/fingerprints/{sample}.json'), 'r') as f:
        data = json.load(f)

    fingerprints = data[0]
    gl_maps = data[1]
    lg_maps = data[2]

    umis = list(sgs_seqs.keys())
    counts_dict = {u: int(u.split('_')[1]) for u in umis}

    G_all = nx.Graph()
    G_all.add_nodes_from(umis)

    G = nx.Graph()
    G.add_nodes_from(umis)
    counts = [int(u.split('_')[1]) / 10 for u in umis]
    # G_trimmed = nx.Graph()
    # G_trimmed.add_nodes_from(umis)
    # counts_trimmed = [int(u.split('_')[1]) / 10 for u in umis]

    G_trimmed = nx.Graph()
    counts_trimmed = []

    G_final = nx.Graph()
    counts_final = []

    print('Detecting fake UMI based on fingerprints...')
    for u1, u2 in combinations(umis, 2):
        umi1 = u1.split('_')[0]
        umi2 = u2.split('_')[0]
        if umi1 == umi2:
            continue
        if seqlib.qhamming(umi1, umi2, limit=1):
            if check_fingerprint_match(umi1, umi2, fingerprints, gl_maps, lg_maps):
                print(u1, u2)
                G.add_edge(u1, u2)
                G_all.add_edge(u1, u2)

    for c in nx.connected_components(G):
        sorted_umi = sorted(c, key=lambda umi: (int(umi.split('_')[1]), umi), reverse=True)
        lumi = sorted_umi[0]
        if len(sorted_umi) > 1:
            for umi in sorted_umi[1:]:
                counts_dict[lumi] += int(umi.split('_')[1])
        G_trimmed.add_node(lumi)
        counts_trimmed.append(counts_dict[lumi] / 10)

    print('Detecting fake UMI based on edit distance and read count...')
    for u1, u2 in combinations(G_trimmed.nodes, 2):
        umi1 = u1.split('_')[0]
        umi2 = u2.split('_')[0]
        if seqlib.qhamming(umi1, umi2, limit=1):
            if int(u1.split('_')[1]) > 2*int(u2.split('_')[1]) or int(u1.split('_')[1]) < 0.5*int(u2.split('_')[1]):
                print(umi1, umi2, int(u1.split('_')[1]), int(u2.split('_')[1]))
                G_trimmed.add_edge(u1, u2)
                G_all.add_edge(u1, u2)

    for c in nx.connected_components(G_trimmed):
        sorted_umi = sorted(c, key=lambda umi: (int(umi.split('_')[1]), umi), reverse=True)
        lumi = sorted_umi[0]
        if len(sorted_umi) > 1:
            for umi in sorted_umi[1:]:
                counts_dict[lumi] += int(umi.split('_')[1])
        G_final.add_node(lumi)
        counts_final.append(counts_dict[lumi] / 10)

    pos = nx.spring_layout(G_all)

    plt.figure(figsize=(25, 10))
    plt.subplot(131)
    nx.draw(G, pos=pos, node_size=counts, linewidths=2, alpha=0.6, with_labels=False)
    plt.title('Fingerprint pairs')

    plt.subplot(132)
    nx.draw(G_trimmed, pos=pos, node_size=counts_trimmed, linewidths=2, alpha=0.6, with_labels=False)
    plt.title('Read count pairs')

    plt.subplot(133)
    nx.draw(G_final, pos=pos, node_size=counts_final, linewidths=2, alpha=0.6, with_labels=False)
    plt.title('Final UMIs')

    plt.savefig(os.path.join(project, f'sgs-prelim/fake-umi/{sample}.png'), format='png', dpi=300)

    with open(os.path.join(project, f'sgs-prelim/fake-umi/{sample}.fasta'), 'w') as f:
        for umi_label in G_final:
            umi = umi_label.split('_')[0]
            final_count = counts_dict[umi_label]
            name = f'{sample}_{umi}_{final_count}'
            f.write(f'>{name}\n{sgs_seqs[umi_label]}\n')