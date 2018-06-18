import core.thread;

import std.algorithm;
import std.array;
import std.container.slist;
import std.conv;
import std.getopt;
import std.json;
import std.parallelism;
import std.path;
import std.random: randomSample;
import std.range;
import std.stdio;
import std.uuid;
import std.zlib: compress;

import mir.random;
import mir.random.engine;
import mir.random.variable;
import mir.math.common: fabs;
import mir.random.algorithm: range, RandomRange, sample;

import vibe.d;


alias RandomRange!(threadLocal!Random, UniformVariable!double) UniformRange;
alias RandomRange!(threadLocal!Random, NormalVariable!double) NormalRange;


const size_t N_TRIALS = 20;
void main (string[] args) {
    
    size_t N_ITER = 1e6.to!size_t;
    size_t nTrials = 10;
    double baseRate = 0.1;
    double fprMutationRate = 0.8;
    double fprMutationMagnitude = 0.05;
    AwardPolicy policy = AwardPolicy.PUBLICATIONS;

    auto helpInformation = getopt(
        args,
        std.getopt.config.passThrough,
        "nTrials", "Number of trials to run", 
            &nTrials,
        "baseRate", "Base rate of true hypotheses", 
            &baseRate,
        "fprMutationRate", "How often the false positive rate mutates", 
            &fprMutationRate,
        "fprMutationMagnitude", "Std. dev. of the false positive mutations",
            &fprMutationMagnitude,
        "policy", "One of: RANDOM, PUBLICATIONS, FPR", 
            &policy
    );

    if (helpInformation.helpWanted)
    {
        defaultGetoptPrinter(
            "SCIMOD\n./scimod-agency WRITE_DIR <OPTIONS>\nOptions:",
            helpInformation.options
        );
        return;
    }
    string writeDir = args[1];

    /****** RUN MANY IN PARALLEL ******/
    defaultPoolThreads(3);

    TrialsData data;
    data.metadata.syncEvery = SYNC_EVERY;
    data.metadata.policy = policy.to!string;
    data.metadata.parameters = [
        "baseRate": baseRate,
        "fprMutationRate": fprMutationRate,
        "fprMutationMagnitude": fprMutationMagnitude,
        "nTrials": nTrials,
        "nIterations": N_ITER
    ];

    data.funds.length = nTrials;
    data.falsePositiveRate.length = nTrials;
    data.nPublications.length = nTrials;

    /* data.funds[0].length = N_ITER / SYNC_EVERY; */
    /* data.falsePositiveRate[0].length = N_ITER / SYNC_EVERY; */
    /* data.nPublications[0].length = N_ITER / SYNC_EVERY; */

    foreach (trialIdx; parallel(nTrials.iota))
    {
        writefln("running %d of %d", trialIdx + 1, nTrials);

        TimeseriesData thisTrialData = simulation(
            policy, baseRate, fprMutationRate, fprMutationMagnitude,
            trialIdx
        );

        data.funds[trialIdx] = thisTrialData.funds;
        data.falsePositiveRate[trialIdx] = thisTrialData.falsePositiveRate;
        data.nPublications[trialIdx] = thisTrialData.nPublications;
    }

    writeDir
        .buildPath(randomUUID().to!string ~ ".json")
        .File("w").write(
            data.serializeToJsonString
        );


    /****** RUN ONE ******/
    /* simulation("test.json"); */
}



const size_t N_PI = 100;
size_t N_ITER = 1e6.to!size_t;
size_t SYNC_EVERY = 2000;
TimeseriesData simulation(AwardPolicy policy,
                double baseRate, double fprMutationRate, 
                double fprMutationMagnitude,
                size_t trialIdx) 
{
    PI[] pis; 

    auto fprMutationAmountRange = 
        normalVar(0.0, fprMutationMagnitude).range;

    auto mutateFprNowRange = uniformVar(0.0, 1.0).range;

    foreach (_; 0..N_PI)
        pis ~= new PI();

    TimeseriesData data;
    data.funds.length = N_ITER / SYNC_EVERY;
    data.falsePositiveRate.length = N_ITER / SYNC_EVERY;
    data.nPublications.length = N_ITER / SYNC_EVERY;

    foreach (iter; 0..N_ITER)
    {
        // PI's try to do research and publish it.
        pis.doScience;
             
        // Agency reviews "grant applications".
        pis.applyForGrants(policy);

        // Select a PI to reproduce and a PI to die.
        pis.evolve(mutateFprNowRange, fprMutationAmountRange, fprMutationRate);
        mutateFprNowRange.popFront();
        fprMutationAmountRange.popFront();

        if (iter % SYNC_EVERY == 0) 
        {
            size_t syncIdx = iter / SYNC_EVERY;

            data.funds[syncIdx].length = N_PI;
            data.falsePositiveRate[syncIdx].length = N_PI;
            data.nPublications[syncIdx].length = N_PI;

            data.funds[syncIdx][] = pis.map!"a.funds".array;
            data.falsePositiveRate[syncIdx][] = 
                pis.map!"a.falsePositiveRate".array;
            data.nPublications[syncIdx][] = pis.map!"a.publications".array;
        }
    }

    return data;
}


private void doScience(PI[] pis) 
{
    foreach (pi; pis)
        pi.doScience();
}


/****************** PRINCIPAL INVESTIGATOR *******************/
const static double SCIENCE_COST = 1.0;
const static double INIT_FUNDS = 100.0;
const double INIT_POWER = 0.5;
const double BASE_RATE = 0.5;
const static double INIT_FALSE_POS_RATE = 0.1;
class PI {
    double funds = INIT_FUNDS;
    double power = INIT_POWER;
    double falsePositiveRate = INIT_FALSE_POS_RATE;
    size_t publications = 0;
    size_t age = 0;
    UniformRange uniformRange;
    
    this() 
    {
        this.uniformRange = uniformVar(0.0, 1.0);
    }

    this(UniformRange uniformRange) 
    {
        this.uniformRange = uniformRange;
    }

    public:
    void doScience()
    {
        if (this.funds > SCIENCE_COST) 
        {
            if (positiveResult()) 
            {
                this.publications += 1;
                this.funds -= SCIENCE_COST;
            }
        }
        ++age;
    }

    void reset()
    {
        this.funds = INIT_FUNDS;
        this.publications = 0;
        this.age = 0;
    }

    private:
        bool positiveResult() {
            return this.uniformRange.front < this.detectionRate;
        }
        // Rate of all detections: fn of base rate, power, and false pos rate.
        @property const double detectionRate() {
            return (BASE_RATE * this.power) + 
                   ((1 - BASE_RATE) * this.falsePositiveRate);
        }
}
unittest {
    PI pi = new PI();
    assert (pi.detectionRate == 0.25 + 0.3);
}


/********* GRANT APPLICATIONS **********/
enum AwardPolicy {RANDOM, PUBLICATIONS}; 


const static size_t N_APPLICANTS = 10;
private void applyForGrants(PI[] pis, AwardPolicy policy=AwardPolicy.RANDOM) 
{
    auto applicants = randomSample(pis, N_APPLICANTS);
    final switch (policy) 
    {
        case AwardPolicy.RANDOM: 

            applicants
                .front
                .addFunds();
            break;

        case AwardPolicy.PUBLICATIONS:

            applicants
                .maxElement!"a.publications"
                .addFunds();
            break;
    }
}


const double AWARD_AMOUNT = 100.0;
private void addFunds(PI pi, double amount=AWARD_AMOUNT)
{
    pi.funds += amount;
}



/********* EVOLUTION **********/
const size_t N_TO_DIE = 10;
const size_t N_TO_REPRODUCE = 10;
private void evolve(
        PI[] pis, 
        UniformRange mutateFprNowRange,
        NormalRange fprMutationAmountRange,
        double fprMutationRate
    )
{
    // "Kill" oldest PI.
    auto resetPI = pis.randomSample(N_TO_DIE).maxElement!"a.age";
    resetPI.reset();
    
    // Reproduce the one with most funding from ten chosen at random.
    auto piToReproduce = pis.randomSample(N_TO_REPRODUCE)
                            .maxElement!"a.publications";

    piToReproduce.reproduce(
        resetPI, mutateFprNowRange, fprMutationAmountRange, fprMutationRate
    );
}


private void reproduce(
        PI piToReproduce, 
        PI resetPI, 
        UniformRange mutateFprNowRange,
        NormalRange fprMutationAmountRange,
        double fprMutationRate
    )
{
    double mutateAmt;
    bool mutate = mutateFprNowRange.front < fprMutationRate;
    
    if (mutate) 
    {
        // Get mutation amount and generate a new random FPR mutation value.
        mutateAmt = fprMutationAmountRange.front;

        // Set FPR of the reset PI as the reproducing PI's FPR plus mutation.
        resetPI.falsePositiveRate = 
            piToReproduce.falsePositiveRate + mutateAmt;
        if (resetPI.falsePositiveRate > 1.0)
        {
            resetPI.falsePositiveRate = 1.0;
        }
        else if (resetPI.falsePositiveRate < 0.0)
        {
            resetPI.falsePositiveRate = 0.0;
        }
    }
}

struct TrialsData
{
    Metadata metadata;
    double[][][] funds;
    double[][][] falsePositiveRate;
    size_t[][][] nPublications;
}


struct TimeseriesData
{
    double[][] funds;
    double[][] falsePositiveRate;
    size_t[][] nPublications;
}


struct Metadata
{
    double baseRate = 0.1;
    size_t syncEvery;
    string policy;
    double[string] parameters;
}
