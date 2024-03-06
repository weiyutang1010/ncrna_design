import numpy as np
import matplotlib.pyplot as plt

with open('./samples_n26.txt', 'r') as f:
    samples = f.read().split('\n')
    samples = [sample for sample in samples if len(sample) > 0]

uniq_seq = set()
num_samples = [(x + 1) * 1000 for x in range(2000)]
unique_seqs = []
prev_num = 0
for num_sample in num_samples:
    for x in samples[prev_num:num_sample]:
        uniq_seq.add(x)
    unique_seqs.append(len(uniq_seq))
    prev_num = num_sample

plt.plot(num_samples, unique_seqs)
plt.xlabel("Number of Samples")
plt.ylabel("Number of Unique Samples")
plt.title(f"Puzzle: ..((((((((.....)).)))))).., n = 26, uniform init, k=1000")
plt.savefig(f'./unique_seq_3.png', format="png", bbox_inches="tight")
