# Sampling-based Continuous Optimization with Coupled Variables for RNA Design
Wei Yu Tang, Ning Dai, Tianshuo Zhou, David H. Mathews, and Liang Huang

## Compile
```
make
```
Compiler version: g++ (Spack GCC) 8.3.0

## Conda Environment (for python scripts)
`conda env create --name ncrna_design --file=environments.yml`

## Scripts
`./main` - run rna design algorithm.

`python analysis.py` - draw learning curve and re-evaluate samples using viennaRNA 2.0

`python excel.py` - compile results and give overall statistics

`./run.sh` - run `./main` in batch

`./graph.sh` - run `python analysis.py` in batch

## Examples
Example commands for running rna design for the shortest 18 structures (up to 50 nucleotides) in Eterna100
```
./run.sh eterna/eterna_n50 example
./graph.sh example
python excel.py eterna_n50 example
```

## ./main
### Example command
```
echo "(((...)))" | ./main [args] > ./result.txt
```

### Arguments

Objective functions: $p(\mathbf{y}^\star \mid \mathbf{x})$ - "prob", $\text{NED}(\mathbf{x}, \mathbf{y}^\star)$ - "ned", $d(\text{MFE}(\mathbf{x}), \mathbf{y}^\star)$ - "dist", $\Delta\Delta G(\mathbf{x}, \mathbf{y}^\star)$ - "ddg". (default: "prob")
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

boxplot: print out the objective of all samples, use to generate boxplot in the learning curve (default: False)
```
--boxplot
```

## ./run.sh
Runs ./main in batch and store results in `./result/[folder]`. Need to manually change the script to add arguments
```
./run.sh <data file> <result folder name>
```

## analysis.py

### Example Command
```
python analysis.py --folder "example" --file "1.txt"
```
Parse the result at `./results/{folder}/{file}` and generate the graph at `./graphs/{folder}/{file}.pdf`

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

## ./graph.sh
Runs `python analysis.py` in batch. Graphs are stored in `./graphs/[folder]/` and best sequences found for each metric are stored in `./analysis/[folder]/`
```
./graph.sh [result folder] <optional: data file>
```

## excel.py
### Example Command
```
python excel.py eterna_n50 example_1 example_2
```
Compile the results from `./analysis/[folder]` and report the average statistics. Can be used to find the best solution out of multiple runs.
```
python excel.py [data file] [result folder 1] <optional: result folder 2, 3, ...>
```