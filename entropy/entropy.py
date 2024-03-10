import numpy as np
import matplotlib.pyplot as plt

with open("../data/eterna/eterna_n50.txt", 'r') as f:
    data = f.read().split('\n')
    data = [(int(line.split(' ')[0]), line.split(' ')[1]) for line in data if len(line) > 0]

for p_id, struct in data:
    with open(f"../results/sampling_uniform_ex/{p_id}.txt", 'r') as f:
        lines = f.read().split('\n')
        n = len(struct)
        num_var = struct.count('.') + struct.count('(')

        steps = [i for i in range(3500)]
        p_entropy = []
        for idx, line in enumerate(lines):
            if line == 'Distribution':
                entropy = []
                for probs in lines[idx+1:idx+1+num_var]:
                    probs = probs.split(': ')[1]
                    probs = np.array([float(prob) for i, prob in enumerate(probs.split(' ')) if i % 2 == 1])
                    temp = 0.
                    for prob in probs:
                        if prob > 0.:
                            temp += prob * np.log2(prob)
                    entropy.append(temp * -1)
                p_entropy.append(np.mean(entropy))

        plt.title("Avg Positional Entropy vs. Steps (Uniform Initialization)")
        plt.xlabel("Steps")
        plt.ylabel("Average Positional Entropy")
        
        linestyle = 'solid'
        # if len(struct) > 31:
        #     linestyle = 'dashed'
        # if len(struct) > 36:
        #     linestyle = 'dotted'

        plt.plot(steps, p_entropy, linestyle=linestyle, label=f'p{p_id}')
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        
plt.savefig(f'./entropy.png', format='png', bbox_inches='tight')