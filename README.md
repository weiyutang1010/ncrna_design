# SamplingDesign: RNA Design via Continuous Optimization with Coupled Variables and Monte-Carlo Sampling

This repository contains the source code for the SamplingDesign project.

Wei Yu Tang, Ning Dai, Tianshuo Zhou, David H. Mathews, and Liang Huang*

\* corresponding author

For questions, please contact the corresponding author at liang.huang.sh@gmail.com.

## To Compile
Compiler version: g++ (Spack GCC) 8.3.0
```
make
```


## Python Dependencies
`python` (3.8.20), `numpy` (1.24.4), `matplotlib` (3.7.5), `viennarna` (2.6.4)
```
pip install -r requirements.txt
```

## Example Script
Run SamplingDesign for the shortest five structures (up to 30 nucleotides) in Eterna100 with 200 steps (takes ~30 seconds).
```
./run.sh example ./data/example.txt
```

The results will be saved in `./results/example/`. The script then parses the result file to generate learning curves in `./graphs/example/` and output the best solution (based on each metric) into `./analysis/example/`.

## To run SamplingDesign
### Command
```
echo "[target structure]" | ./main [args] > ./result.txt
```
### Example
```
echo "(((...)))" | ./main --steps 50 --verbose > ./result.txt
```

### Arguments

Objective functions: "prob" - Boltzmann probability, "ned" - normalized ensemble defect, "dist" - structural distance, "ddg" - free energy gap. (default: "prob")
```
--obj [prob/ned/dist/ddg]
```

Initializations: uniform or targeted (default: targeted)
```
--init [uniform/targeted]
```

eps: For $\epsilon$-targeted initialization. To use targeted initialization, set $\epsilon$ = 1.0.  (default: 0.75)
```
--eps [a value between 0 and 1]
```

projection: use the direct parameterization (projected gradient descent) instead of the softmax parameterization. (default: false)
```
--projection
```

no_adam: turn off adam optimizer (default: false)
```
--no_adam
```

beta_1, beta_2: the first and second moment decay rates for adam optimizer (default: 0.9, 0.999)
```
--beta_1 [value] --beta_2 [value]
```

nesterov: use Nesterov's accelerated gradient descent for projected gradient descent (default: false)
```
--nesterov
```

initial_lr: set initial learning rate (default: 0.01)
```
--initial_lr [value]
```

<!-- lr_decay, lr_decay_rate: use learning rate decay (default: false, 0.96)
```
--lr_decay --lr_decay_rate [value]
```

adaptive_lr, k_ma_lr: use adaptive learning rate decay (default: false, 10)
```
--adaptive_lr --k_ma_lr [value]
``` -->

num_steps: set max number of steps (default: 2000)
```
--num_steps [value]
```

adaptive_step: use early stopping condition (default: true, 50)
```
--adaptive_step
```

k_ma: k moving average
```
--k_ma [value]
```

beamsize: set beamsize and sharpturn (default: 250, False)
```
--beamsize [value]
```

is_lazy: use lazyOutside when evaluating NED (default: False)
```
--is_lazy
```

sample_size: number of samples used per step (default: 2500)
```
--sample_size [value] 
```

best_k: print out best k samples (in terms of objective function) at each step (default: 1)
```
--best_k [value]
```

<!-- no_mismatch, no_trimismatch: turn off coupling variables for mismatch and trimismatch (default: False, False)
```
--no_mismatch, --no_trimismatch
``` -->

verbose: print out the logit, the distribution and the gradient at each step (default: False)
```
--verbose
```
<!-- 
seed: random seed (default: 42)
```
--seed [value]
``` -->

num_threads: max number of threads used by openMP, if 0 then use the default number (default: 0)
```
--num_threads [value]
```

boxplot: print out the objective of all samples at each step (for generating the boxplot in the learning curve) (default: False)
```
--boxplot
```



## analysis.py
`analysis.py` performs two main tasks:
1. Re-evaluates the best solution at each step using the ViennaRNA package (version 2.6.4), and saves the best solution for each metric to ./analysis/{folder}/.
2. Generates a learning curve and saves it to ./graphs/{folder}/.

### Example Command
```
python analysis.py --folder "example" --file "1.txt"
```
Input: Parses the results from ./results/{folder}/{file}.
Output:

- Best solutions are saved to ./analysis/{folder}/{file}.txt

- The learning curve is saved as ./graphs/{folder}/{file}.pdf

### Arguments
folder, file: parse the file at `./results/{folder}/{file}`
```
--folder "example" --file "8.txt"
```

max_workers: set the max number of workers (threads) (default: None)
```
--max_workers [value]
```
<!-- 
no_eval_seq: produce learning curve only (default: False)
```
--no_eval_seq
```

no_graph: re-evaluate all sequences only (default: False)
```
--no_graph
``` -->

## Reproduction Instructions
Reproduce the results in the paper (uniform and $\epsilon$-targeted initializations). Note that it can take longer than 2 weeks to run the entire script.

### Command
```
./run_all.sh ./data/eterna100.txt
```
