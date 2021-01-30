Here are some example of runs:

Unfolded data, all models:

```bash
grapes -in Microtus.unfolded.dofe -out Microtus.unfolded -model all
```

Folded data, all models:

```bash
grapes -in Microtus.folded.dofe -out Microtus.folded -model all
```


Running all models on folded data, from an unfolded input SFS:

```bash
grapes -in Microtus.unfolded.dofe -out Microtus.folded2 -model all -fold
```

Running only a subset of models:

```bash
grapes -in Microtus.unfolded.dofe -out Microtus.unfolded2 -model DisplGamma,GammaExpo,ScaledBeta
```



