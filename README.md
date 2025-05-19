# SamplingDesign: RNA Design via Continuous Optimization with Coupled Variables and Monte-Carlo Sampling
Wei Yu Tang, Ning Dai, Tianshuo Zhou, David H. Mathews, and Liang Huang

## Compile
```
make
```
Compiler version: g++ (Spack GCC) 8.3.0

## Conda Environment (for python scripts)
`conda env create --name ncrna_design --file=environments.yml`

## Script
Run SamplingDesign for the shortest five structures (up to 30 nucleotides) in Eterna100 (100 steps).
```
./run.sh example ./data/example.txt
```

```
./run.sh <result folder> <data file path>
```
The results are stored in `./results/`. The script also output learning curves in `./graphs/` and solution in `./analysis/`

## ./main
### Command
```
echo "[target structure]" | ./main [args] > ./result.txt
```
Example:
```
echo "(((...)))" | ./main --steps 50 --verbose > ./result.txt
```

### Arguments

Objective functions: "prob" - Boltzmann probability, "ned" - normalized ensemble defect, "dist" - structural distance, "ddg" - free energy gap. (default: "prob")
```
--obj [prob/ned/dist/ddg]
```

Initializations: uniform/targeted/eps-targeted (default: targeted)
```
--init [uniform/targeted]
```

eps: used for eps-targeted initialization. for targeted, set eps = 1.0.  (default: 0.75)
```
--eps [value between 0 and 1]
```

projection: use projection instead of softmax parameterization. (default: false)
```
--projection
```

no_adam: turn off adam optimizer (default: false)
```
--no_adam
```

nesterov: use Nesterov's accelerated gradient descent (default: false)
```
--nesterov
```

beta_1, beta_2: first and second moment decay rates for adam optimizer (default: 0.9, 0.999)
```
--beta_1 [value] --beta_2 [value]
```

initial_lr: set initial learning rate (default: 0.01)
```
--initial_lr [value]
```

lr_decay, lr_decay_rate: use learning rate decay (default: false, 0.96)
```
--lr_decay --lr_decay_rate [value]
```

adaptive_lr, k_ma_lr: use adaptive learning rate decay (default: false, 10)
```
--adaptive_lr --k_ma_lr [value]
```

num_steps: set max number of steps (default: 2000)
```
--num_steps [value]
```

adaptive_step, k_ma: use adaptive stopping condition (default: false, 50)
```
--adaptive_step --k_ma [value]
```

beamsize, sharpturn: set beamsize and sharpturn (default: 250, False)
```
--beamsize [value] --sharpturn
```

is_lazy: use lazyOutside when evaluating NED (default: False)
```
--is_lazy
```

sample_size, best_k: number of samples and print out top k sample at each step (default: 2500, 1)
```
--sample_size [value] --best_k [value]
```

no_mismatch, no_trimismatch: turn off coupling variables for mismatch and trimismatch (default: False, False)
```
--no_mismatch, --no_trimismatch
```

verbose: print out logit, distribution and gradient at each step (default: False)
```
--verbose
```

seed: random seed (default: 42)
```
--seed [value]
```

num_threads: max number of threads used by openMP, if 0 then use default number (default: 0)
```
--num_threads [value]
```

boxplot: print out the objective of all samples at each step (for generating the boxplot in the learning curve) (default: False)
```
--boxplot
```



## analysis.py

### Example Command
```
python analysis.py --folder "example" --file "1.txt"
```
Input: takes the results from `./results/{folder}/{file}`. Output: list best solutions at `./analysis/{folder}/{file}.txt` and generate the learning curve at `./graphs/{folder}/{file}.pdf`

### Arguments
folder, file: parse the file at `./results/{folder}/{file}`
```
--folder "example" --file "8.txt"
```

max_workers: set the max number of workers (default: None)
```
--max_workers [value]
```

no_eval_seq: produce learning curve only (default: False)
```
--no_eval_seq
```

no_graph: reeval all sequences only (default: False)
```
--no_graph
```
