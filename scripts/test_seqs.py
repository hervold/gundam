#!/usr/bin/env python3

from argparse import ArgumentParser
import random
import sys

SEQLEN = 500
ERRRATE = 0.02

def generate_seqs(args):
    prob = 0.05
    motifs = []
    for idx in range(args.motif_ct):
        m = ''.join(random.choice('ATGC') for _ in range(random.randint(10,20)))
        motifs.append((m, prob))
        prob = prob / 2.0
        print('>motif_{}'.format(idx))
        print(m)

    for idx in range(args.seq_ct):
        seq = ''.join(random.choice('ATGC') for _ in range(args.seqlen))
        u = random.random()
        #sys.stderr.write('-- {}\n'.format(u))
        for m, p in motifs:
            if u <= p:
                loc = random.randint(0, args.seqlen - len(m) - 1)
                m2 = ''.join(random.choice('ATGC') if random.random() <= ERRRATE else m[i] for i in range(len(m)))
                seq = seq[:loc] + m2 + seq[loc+len(m2):]
        print('>x{}'.format(idx))
        print(seq)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('seq_ct', type=int, help='number of sequences to generate')
    parser.add_argument('motif_ct', type=int, help='number of distinct motifs')
    parser.add_argument('--seqlen', type=int, default=SEQLEN, help='sequence length (default {})'.format(SEQLEN))

    args = parser.parse_args()
    generate_seqs(args)
