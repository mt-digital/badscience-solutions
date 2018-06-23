import core.thread;

import std.algorithm;
import std.array;
import std.container.slist;
import std.conv;
import std.getopt;
import std.json;
import std.math: abs, approxEqual;
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

import vibe.d: serializeToJsonString;


alias RandomRange!(threadLocal!Random, UniformVariable!double) UniformRange;
alias RandomRange!(threadLocal!Random, NormalVariable!double) NormalRange;


const size_t N_TRIALS = 20;
int main (string[] args) {

    size_t N_ITER = 1e6.to!size_t;
    size_t nTrials = 10;
    double awardAmount = 50.0;
    double seedFunding = 10.0;
    double baseRate = 0.1;
    double fprMutationRate = 0.25;
    double fprMutationMagnitude = 0.01;
    double initialFalsePositiveRate = 0.05;
    double publishNegativeResultRate = 0.0;
    double falsePostivePaperDiscoveryRate = 0.0;
    double grantApplicationCost = 0.0;
    AwardPolicy policy = AwardPolicy.FPR;

    auto helpInformation = getopt(
        args,
        std.getopt.config.passThrough,
        "nTrials", "Number of trials to run (default 10)", 
            &nTrials,
        "baseRate", "Base rate of true hypotheses (default 0.1)", 
            &baseRate,
        "awardAmount", "Amount given to grant-winning lab in a timestep (default 50)", 
            &awardAmount,
        "initialFalsePositiveRate", "False positive rate of all PIs at t=0 (default 0.05)", 
            &initialFalsePositiveRate,
        "fprMutationRate", "How often the false positive rate mutates (default 0.25)", 
            &fprMutationRate,
        "fprMutationMagnitude", "Std. dev. of the false positive mutations (default 0.01)",
            &publishNegativeResultRate,
        "publishNegativeResultRate", "Rate that negative results are published (default 0.0)",
            &fprMutationMagnitude,
        "policy", "One of: RANDOM, PUBLICATIONS, FPR (default PUBLICATIONS)", 
            &policy
    );

    if (helpInformation.helpWanted)
    {
        defaultGetoptPrinter(
            "SCIMOD\n./scimod-agency WRITE_DIR <OPTIONS>\nOptions:",
            helpInformation.options
        );
        return 0;
    }

    if (args.length == 1) {
        writeln("Write directory not provided!");
        return 1;
    }
    string writeDir = args[1];

    /****** RUN MANY IN PARALLEL ******/
    defaultPoolThreads(totalCPUs);
    TrialsData data;
    data.metadata.syncEvery = SYNC_EVERY;
    data.metadata.policy = policy.to!string;
    data.metadata.parameters = [
        "baseRate": baseRate,
        "fprMutationRate": fprMutationRate,
        "fprMutationMagnitude": fprMutationMagnitude,
        "awardAmount": awardAmount,
        "nTrials": nTrials,
        "nIterations": N_ITER,
        "publishNegativeResultRate": publishNegativeResultRate
    ];

    data.funds.length = nTrials;
    data.falsePositiveRate.length = nTrials;
    data.nPublications.length = nTrials;

    foreach (trialIdx; parallel(nTrials.iota))
    {
        writefln("running %d of %d", trialIdx + 1, nTrials);

        TimeseriesData thisTrialData = simulation(
            policy, awardAmount, baseRate, initialFalsePositiveRate,
            fprMutationRate, fprMutationMagnitude, publishNegativeResultRate
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

    return 0;
}



const size_t N_PI = 100;
size_t N_ITER = 1e6.to!size_t;
size_t SYNC_EVERY = 2000;
TimeseriesData simulation(AwardPolicy policy,
                double awardAmount, double baseRate, 
                double initialFalsePositiveRate, double fprMutationRate, 
                double fprMutationMagnitude, double publishNegativeResultRate) 
{
    PI[] pis; 

    auto fprMutationAmountRange = 
        normalVar(0.0, fprMutationMagnitude).range;

    auto mutateFprNowRange = uniformVar(0.0, 1.0).range;

    foreach (i; 0..N_PI)
    {
        pis ~= new PI();
        pis[i].funds = awardAmount;
        pis[i].falsePositiveRate = initialFalsePositiveRate;
        pis[i].publishNegativeResultRate = publishNegativeResultRate;
    }

    TimeseriesData data;
    data.funds.length = N_ITER / SYNC_EVERY;
    data.falsePositiveRate.length = N_ITER / SYNC_EVERY;
    data.nPublications.length = N_ITER / SYNC_EVERY;

    foreach (iter; 0..N_ITER)
    {
        // PI's try to do research and publish it.
        pis.doScience;
             
        // Agency reviews "grant applications".
        pis.applyForGrants(awardAmount, policy);

        // Select a PI to reproduce and a PI to die.
        pis.evolve(mutateFprNowRange, fprMutationAmountRange, 
                   awardAmount, fprMutationRate);
        mutateFprNowRange.popFront();
        fprMutationAmountRange.popFront();

        if (iter % SYNC_EVERY == 0) 
        {
            size_t syncIdx = iter / SYNC_EVERY;

            data.funds[syncIdx] = 
                pis.map!"a.funds".array.mean;

            data.falsePositiveRate[syncIdx] = 
                pis.map!"a.falsePositiveRate".array.mean;

            data.nPublications[syncIdx] = 
                pis.map!"a.publications".array.mean;
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
const static double INIT_FUNDS = 10.0;
const double INIT_POWER = 0.8;
const double BASE_RATE = 0.1;
const static double INIT_FALSE_POS_RATE = 0.05;
class PI {
    double funds = INIT_FUNDS;
    double power = INIT_POWER;
    double falsePositiveRate = INIT_FALSE_POS_RATE;
    double publishNegativeResultRate = 0.0;
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
            if (this.funds >= SCIENCE_COST) 
            {
                this.funds -= SCIENCE_COST;
                if (positiveResult()) 
                {
                    // TODO: if it's a false positive, have it be dected with r_d
                    // and do not publish
                    this.publications += 1;
                }
                else if (publishNegativeResult()) 
                { 
                    this.publications += 1;
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
        bool positiveResult() 
        {
            bool ret = this.uniformRange.front < this.detectionRate;
            this.uniformRange.popFront();
            return ret;
        }

        bool publishNegativeResult() 
        {
            bool ret = this.uniformRange.front < this.publishNegativeResultRate;
            this.uniformRange.popFront();
            return ret;
        }
        // Rate of all detections: fn of base rate, power, and false pos rate.
        @property const double detectionRate() {
            return (BASE_RATE * this.power) + 
                   ((1 - BASE_RATE) * this.falsePositiveRate);
        }
}
unittest {
    PI pi = new PI();
    assert(pi.detectionRate.approxEqual(0.125));

    // Set up PI with enough funds to do 100 rounds of research.
    pi.funds = 100;
    foreach (_; 0..100)
    {
        pi.doScience();
    }
    // Statistically this should be true almost all the time.
    assert(pi.publications > 5, pi.publications.to!string);
    assert(pi.age == 100);
    assert(pi.funds == 0, pi.funds.to!string);

    pi.reset();
    assert(pi.publications == 0);
    assert(pi.age == 0);
    assert(pi.funds == INIT_FUNDS);

    // Now test publishing negative results.
    // Make it impossible for the PI to publish positive results and count
    // negative results.
    pi.power = 0.0;
    pi.falsePositiveRate = 0.0;
    pi.publishNegativeResultRate = 0.5;
    pi.funds = 1000;
    foreach (_; 1000.iota)
    {
        pi.doScience(); 
    }
    assert(abs(pi.publications - 500.0) < 50.0, pi.publications.to!string);
}


/********* GRANT APPLICATIONS **********/
enum AwardPolicy {RANDOM, PUBLICATIONS, FPR}; 


const static size_t N_APPLICANTS = 10;
private void applyForGrants(PI[] pis, double awardAmount, 
                            AwardPolicy policy=AwardPolicy.RANDOM,
                            double epsilon=0.1) 
{
    auto applicants = randomSample(pis, N_APPLICANTS);
    final switch (policy) 
    {
        case AwardPolicy.RANDOM: 

            applicants
                .front
                .addFunds(awardAmount);
            break;

        case AwardPolicy.PUBLICATIONS:

            applicants
                .maxElement!"a.publications"
                .addFunds(awardAmount);
            break;

        case AwardPolicy.FPR:

            applicants
                .minElement!"a.falsePositiveRate"
                .addFunds(awardAmount);
            break;
    }
}
unittest 
{
    PI[] pis;
    foreach (i; 0..10)
    {
        pis ~= new PI();
        pis[i].publications = i;
        pis[i].falsePositiveRate = (i + 1) * .05;
    }

    pis.applyForGrants(10.0, AwardPolicy.PUBLICATIONS);
    pis.applyForGrants(10.0, AwardPolicy.FPR);
    double fprMinFunds = pis[0].funds;
    assert(fprMinFunds == INIT_FUNDS + 10.0);

    double pubMaxFunds = pis[$-1].funds;
    assert(pubMaxFunds == INIT_FUNDS + 10.0);
}


// TODO make PI method.
private void addFunds(PI pi, double amount)
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
        double initialFunds,
        double fprMutationRate
    )
{
    // "Kill" oldest PI.
    auto resetPI = pis.randomSample(N_TO_DIE).maxElement!"a.age";
    resetPI.reset();
    resetPI.funds = initialFunds;
    
    // Reproduce the one with most funding from ten chosen at random.
    auto piToReproduce = pis.randomSample(N_TO_REPRODUCE)
                            .maxElement!"a.publications";

    piToReproduce.reproduce(
        resetPI, mutateFprNowRange, fprMutationAmountRange, fprMutationRate
    );
}
unittest
{
    PI[] pis;
    foreach (i; 0..10)
    {
        pis ~= new PI();
        pis[i].age = i;
        pis[i].publications = (10 - i);
        pis[i].falsePositiveRate = 0.1;
    }
    pis[0].falsePositiveRate = 0.314159;

    double fprMutationMagnitude = 0.01;
    auto fprMutationAmountRange = 
        normalVar(0.0, fprMutationMagnitude).range;

    auto mutateFprNowRange = uniformVar(0.0, 1.0).range;
    double awardAmount = 1.0;
    double fprMutationRate = 0.0;
    pis.evolve(
        mutateFprNowRange, fprMutationAmountRange, awardAmount, 
        fprMutationRate
    );
    mutateFprNowRange.popFront();
    fprMutationAmountRange.popFront();

    // The oldest of the ten PIs should die and inherit the FPR of the
    // PI with the most publications.
    assert(pis[$-1].age == 0);
    assert(pis[$-1].falsePositiveRate == 0.314159);

    fprMutationRate = 0.9;
    foreach (pi; pis)
        pi.falsePositiveRate = 0.1;

    foreach (i; 0..5)
    {
        pis.evolve(
            mutateFprNowRange, fprMutationAmountRange, awardAmount, 
            fprMutationRate
        );
        mutateFprNowRange.popFront();
        fprMutationAmountRange.popFront();
    }
    assert(pis.map!"a.falsePositiveRate".sum != 1.0);
}


// TODO make this a PI method.
private void reproduce(
        PI piToReproduce, 
        PI resetPI, 
        UniformRange mutateFprNowRange,
        NormalRange fprMutationAmountRange,
        double fprMutationRate
    )
{
    double mutateAmt = 0.0;
    bool mutate = mutateFprNowRange.front < fprMutationRate;
    
    if (mutate) 
    {
        // Get mutation amount and generate a new random FPR mutation value.
        mutateAmt = fprMutationAmountRange.front;
    }

    // Set FPR of the reset PI as the reproducing PI's FPR plus mutation.
    resetPI.falsePositiveRate = piToReproduce.falsePositiveRate + mutateAmt;

    if (resetPI.falsePositiveRate > 1.0)
    {
        resetPI.falsePositiveRate = 1.0;
    }
    else if (resetPI.falsePositiveRate < 0.0)
    {
        resetPI.falsePositiveRate = 0.0;
    }
}

struct TrialsData
{
    Metadata metadata;
    double[][] funds;
    double[][] falsePositiveRate;
    double[][] nPublications;
}


struct TimeseriesData
{
    double[] funds;
    double[] falsePositiveRate;
    double[] nPublications;
}


struct Metadata
{
    double baseRate = 0.1;
    size_t syncEvery;
    string policy;
    double[string] parameters;
}
