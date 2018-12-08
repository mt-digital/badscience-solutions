# scimod-agency

What can science funding agencies do to improve the quality of the
scientific work they fund? They can change how much they give to award
recipients, or the policy of whom to give awards to. We show that the proper
choice of how much and to whom depends on cultural factors of a scientific
field: how often negative results are published and how good is peer review.

## run the model

### build the code

To run the model first you must get and build the code. To do that you need
to get the D compiler installed, called `dmd`.  On macOS, type `brew install dmd`
into the terminal and press enter. We use
[dub](http://code.dlang.org/getting_started) to build our compiled executable, 
so get dub installed by following the instructions at the link.

When all this is done, fetch the code by cloning this repository,
then `cd scimod-agency` and build the code by running `dub build`. Run the
unit tests using `dub test`. 

### run executable `./scimod-agency`

The `scimod-agency` executable prints its help like so

```
./scimod-agency -h

SCIMOD
./scimod-agency WRITE_DIR <OPTIONS>
Options:
                      --nTrials Number of trials to run (default 10)
                     --baseRate Base rate of true hypotheses (default 0.1)
                  --awardAmount Amount given to grant-winning lab in a timestep (default 50)
     --initialFalsePositiveRate False positive rate of all PIs at t=0 (default 0.05)
              --fprMutationRate How often the false positive rate mutates (default 0.25)
    --publishNegativeResultRate Rate that negative results are published (default 0.0)
         --fprMutationMagnitude Std. dev. of the false positive mutations (default 0.01)
   --falsePositiveDetectionRate Std. dev. of the false positive mutations (default 0.01)
                       --policy One of: RANDOM, PUBLICATIONS, FPR (default PUBLICATIONS)
-h                       --help This help information.
```

#### Build the parameters files

There is a command line option to pass a comma-separated tuple of the four
parameters we varied in our experiments. The four in order are policy, award
amount, negative results publishing rate, and false positive detection rate.
This option is meant to be used in conjunction with the script to make files
with one parameter tuple on each line, `experiment_makeparams.sh`:

```bash
./experiment_makeparams.sh | cat > finaldraft-params.txt
```

1452 parameter combinations are contained in `finaldraft-params.txt`.
Because the MERCED cluster limits the number of jobs allowed in a job array to
1000, we split the `finaldraft-params.txt` in two: 

```bash
tail -n452 finaldraft-params.txt > finaldraft-params-2.txt
```

#### Deploy job arrays to queue

Job arrays make submitting many parameter combinations easy. With the parameter
files set up using the above instructions, we submit the jobs to the cluster
in two steps:

```bash
sbatch --array=1-1000 experiment.sh finaldraft-params.txt
sbatch --array=1-452 experiment.sh finaldraft-params-2.txt
```

It is somewhat sloppy, but for other parameter sensitivity analysis, we just
commented out a block of code in `experiment.sh`

## Data pipeline

Currently there is a process that must be done to convert the directory of
JSONs created by the distributed model runs into a single HDF. It should be
changed to be submitted to the scheduler asap; it can be multithreaded. See
`sandbox/multiprocessing_example.py`.  

To run,

```sh
python json_to_hdf.py path/to/jsons/dir new.hdf
```

Doing this on the cluster, you need to submit it as a job to the queue, like so

```bash
sbatch json_to_hdf_slurm.sh ~/scr/scimod-baseRate0.5/ ~/scr/scimod-baseRate0.5.hdf
```
for example.


The HDF can be read as an `ExperimentData` instance

```python
from experiment_data import ExperimentData
ed = ExperimentData('new.hdf')
```

This object takes advantage of Python indexing ordered as `policy`,
`award_amount`, `pubneg_rate`, and `fpdr`. For example, 

```python
policy = 'FPR'; award_amount = 10; pubneg_rate = '0.10'; fpdr = '0.10'
measures = ed[policy, award_amount, pubneg_rate, fpdr]
assert list(measures.keys()) == [
    'falseDiscoveryRate', 'falsePositiveRate', 'meanFunds', 'meanPublications', 
    'medianFunds', 'medianPublications', 'sumFunds', 'sumPublications'
]
```
