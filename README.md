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

By writing scripts like `run_negres-peer-rev_experiment.sh`, shown below,
we test many parameter combinations.

```bash
for policy in RANDOM PUBLICATIONS; do
    for awardAmount in 1 `seq 5 5 115` ; do
        for publishNegativeResultRate in 0.5 0.9; do
            for falsePositiveDetectionRate in 0.5 0.9; do
                echo \
                    "policy=$policy; awardAmount=$awardAmount;"\
                    "pubNegResRate=$publishNegativeResultRate;" \ 
                    "falsePositiveDetectionRate=$falsePositiveDetectionRate"

                qsub -S /bin/bash -q fast.q -cwd -j y -V -l mem_free=96G -pe \
                    smp 20 -N negres-peer -o negres-peer.log -e negres-peer.err \
                    negres-peer-rev_trials_qsub.sh $policy $awardAmount \
                    $publishNegativeResultRate $falsePositiveDetectionRate

        done
        done
    done
done
```

## Data pipeline

Currently there is a process that must be done to convert the directory of
JSONs created by the distributed model runs into a single HDF. It should be
changed to be submitted to the scheduler asap; it can be multithreaded. See
`sandbox/multiprocessing_example.py`.  

To run,

```sh
python json_to_hdf.py path/to/jsons/dir new.hdf
```
