```python plot.py [-p] [-c] [-l] [--mfe] [-t #] [-o <folder>]```

```
Mode
-p: draw positional defects plots
-c: draw circular bpp plots
-l: draw linear bpp plots
--mfe: draw mfe plots
-t #: use layout # from the viennaRNA package
-o <folder>: save plots to ./plots/<folder> (default folder is "temp")
```

The scripts read from the ```input.txt``` in the format of
```
id1
id2
...
seq1
seq2
...
struct1
struct2
...
```