# ExpectedPartition: RNA Design with Projected Gradient Descent


## To Compile
```
make
```

## Files
src/ExpectedPartition.cpp - Contains code for main, gradient_descent, projection

src/ExpectedPartition.h - Contains code for class definitions and function headers

src/Inside.cpp - Contains code for expected inside partition

src/Outside.cpp - Contains code for expected outside partition (gradient)

src/Eval.cpp - Contains code for calculating Delta_G(x, y) and eval() function

## To Run
ExpectedPartition can be run with:
```
echo STRUCTURE | ./expectedpartition [OPTIONS]
```

OPTIONS:
flags:
```
--obj
```
Objective Function, 0: -log p(y|x), 1: Delta_G(x, y) (default: '0')

```
--init
```
initialization: 0 - stdin, 1 - uniform, 2 - target (default: '2')

```
--lr
```
set learning rate (default: '0.001')

```
--step
```
set number of iterations (default: '1000')

```
--penalty
```
penalty value

```
--eval
```
evaluate p(y|x) given the seq and struct (default: 'false')

## To Test
```
make test
```
or
```
cat data/eterna/short_eterna.txt | ./expectedpartition --init 2 --obj 0 --test --verbose
```