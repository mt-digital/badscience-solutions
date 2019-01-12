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


## Subtler funding strategies

In our first round of experiments, we tested only very simple strategies the
funding agency would use for deciding which PI received a grant: totally
random, given to PI with best methodological integrity in a random sample of ten
as evidenced by lowest false positive rate, 
and given to the PI with the most publications in a random sample of ten. 
In addition to these, we considered the _modified random_ and _mixed_ funding
strategies. In the _modified random_ strategy, a PI is only qualified for
receiving a grant if its false positive rate is not greater than some threshold
value we call _A_ in the paper. In that strategy the grant is awarded at
random to one qualified PI. In the _mixed_ strategy, the grant is awarded to
the PI with the best methodological integrity a fraction $X$ of the time, with
the grant awarded randomly the other $1-X$ of the time. As a reminder,
best methodological integrity means the PI with the minimum false positive rate
among all PIs.

These experiments required us to add an additional parameter, `--policyParam`,
set to zero by default. It is only used if one of the two new policies are
specified, indicated by strings `MODIFIED_RANDOM` and `MIXED`. We also enabled
the `policyParam` to be passed as a fifth comma-separated value in the 
`paramsList` option, e.g. `--paramsList=MIXED,85,0.3,0.3,0.8` says use the
mixed funding strategy, funding per grant of 85, peer review efficacy and
publication rate of negative results of 0.3, and the `policyParam`, $X$, of 0.8.

In our experiments with these we generate .txt lists of comma-separated parameters
used in submitting an array job to the Slurm cluster. Generate the .txt files
and submit the jobs using the following lines of code. Check that the length of
the .txt files is 484. Make an appropriately-named directory in scratch, follow
the examples in `experiment.sh` and put that directory name in the appropriate
place. There is certainly a better way to streamline all this; we or others
will do that later if it makes sense to.

```bash
# Build parameter lists, one for each funding strategy.
./experiment_makeparams_suppPolicies_modran.sh > modran-final-params.txt
./experiment_makeparams_suppPolicies_mixed.sh > mixed-final-params.txt

# Submit jobs to the cluster.
sbatch --array=1-484 experiment.sh modran-final-params.txt
sbatch --array=1-484 experiment.sh mixed-final-params.txt
```

Then when it's time to analyze the results, there is not yet a converter to
HDF. Nonetheless, the process is zippy. `scp` the two directories with 484
JSONs. Let's say you do that and now the two directories are `modran-dir` and
`mixed-dir`. Here's how to load the JSONs from each of these directories into
one single JSON, then pass with the right auxiliary arguments to the heatmap
plotting routine. I just executed this in an IPython shell.
Loading each JSON takes well over a minute on my MacBook Pro. 

```python
from vis import all_supplemental_policy_heatmaps, _make_json_dict

# Pre-process 2x484 JSON: calculate avg final FPR and final FDR across dims.
jd_modran = _make_json_dict('modran-dir')
jd_mixed = _make_json_dict('mixed-dir')

# Peer review efficacy and negative publishing rates, equal in this experiment.
fpdr_npr_rates = np.arange(0.0, 1.01, 0.1)

# Can specify different scales for each strategy. We use the same scale here.
policy_params_dict = {
    'MODIFIED_RANDOM': fpdr_npr_rates, 
    'MIXED': fpdr_npr_rates
}


# Could have also created JSONs using this method, so it returns them back
# unchanged. This will create 16 heatmap figures saved to 
# os.path.expanduser('~/workspace/papers/sciencefunding/Figures/') -- not 
# useful to most. You might have to make some changes to change this behavior
# more easily. 16 = 4 funding per grant values x 2 funding strategies x 2
# measures of research quality (ave FPR and FDR).

jd_modran, jd_mixed = all_supplemental_policy_heatmaps(
    json_dict_modran=jd_modran, json_dict_mixed=jd_mixed,
    fpdr_npr_rates=fpdr_npr_rates, policy_params_dict=policy_params_dict
)
```

## Other supplemental analyses, figures, etc.

We have generated a supplement that demonstrates model convergence
and that our choice of base rate and selection process do not influence 
our results. Below I briefly explain these checks further as I show you how
the checks were done and visualized using the software in this repository.

One colleague who graciously reviewed a preliminary version of the paper claimed
his base rate was 0.5. On our view, a base rate of 0.1 may be inflated for
most fields/researchers. In any case, we thought since 0.5 is on the extreme
end of reasonable, we'll use that. For testing this parameter setting, we
clumsily copy/pasted the single-line block of a bash command in `experiment.sh`, 
commented out the original, and added the option `--baseRate=0.5` to the 
command block. See `experiment.sh` for this command, which itself is now 
commented out.

Another colleague helpfully suggested we use an alternative selection method.
Our original selection method was to select ten PIs at random, then the one
with the most publications of the ten reproduced. This alternative selection
method our colleague called "Wright-Fisher", so we did too. Here we see the
beauty of D in the implementation of Wright-Fisher selection in `source/app.d`:

```d
reproducingIdx = dice(pis.map!"a.publications");
```

So, the PI is selected at random with probability equal to the number of its
publications divided by the sum of publications over all PIs. A short aside 
on what's going on to help with D: `map` is not a 
method of `pis`, which is just an array of `PI` instances, i.e. `PI[]`. No,
instead this is D's universal function call syntax (UFCS), where `pis` is actually the
first (and in this case only) argument to `map`, 
which can be found in the standard library's
`std.algorithm` module. Since `pis` is the only argument, no parentheses are 
needed at the end of this function call in D. If there were more arguments,
they would go in parentheses after the second quote.
The bang, `!`, indicates the start of a template
argument, which can be a string that defines a function, in this case 
`"a.publications"`. Making alternative choices for both UFCS and string
definition of the anonymous function template argument, we get, 

```d
reproducingIdx = dice(map!(anonFuncVar => anonFuncVar.publications)(pis));
```

where `anonFuncVar` is arbitrary. In the string version, one must
use `a` as the first function argument.

Back to the experiment in the paper, we tested Wright-Fisher selection
with a base rate of 0.1 and 0.5. These commented-out calls to `./scimod-agency` 
(the old name of the repo/program) can be found in `experiment.sh`. To use
the Wright-Fisher selection instead of the best of random ten, add the option
`--selectionMethod=WRIGHT_FISHER` to your `./scimod-agency` call.
