#!/usr/bin/env python3
import argparse, numpy as np


def hist_text(vals, bins=20):
vals = np.asarray(vals)
if vals.size == 0:
return "(no reads)\n"
counts, edges = np.histogram(vals, bins=bins)
lines = []
for i,(c,lo,hi) in enumerate(zip(counts, edges[:-1], edges[1:])):
lines.append(f"bin{i+1:02d}\t{int(lo)}-{int(hi)}\t{int(c)}")
return "\n".join(lines) + "\n"


def main():
ap = argparse.ArgumentParser()
ap.add_argument('lengths')
ap.add_argument('--out-all', required=True)
ap.add_argument('--out-mid80', required=True)
ap.add_argument('--out-low80', required=True)
ap.add_argument('--bins', type=int, default=20)
a = ap.parse_args()
with open(a.lengths) as fh:
L = [int(x.strip()) for x in fh if x.strip()]
n = len(L)
with open(a.out_all,'w') as o:
o.write(hist_text(L, a.bins))
# middle 80% (drop 10% head/tail)
Ls = sorted(L)
k = int(round(n*0.1))
mid = Ls[k:n-k] if n>=10 else Ls
with open(a.out_mid80,'w') as o:
o.write(hist_text(mid, a.bins))
# low 80% (exclude top 20%)
k2 = int(round(n*0.2))
low80 = Ls[:n-k2] if n>=5 else Ls
with open(a.out_low80,'w') as o:
o.write(hist_text(low80, a.bins))


if __name__ == '__main__':
main()