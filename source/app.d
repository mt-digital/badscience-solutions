import core.thread;

import std.algorithm;
import std.array;
import std.container.slist;
import std.conv;
import std.exception;
import std.getopt;
import std.json;
import std.math: abs, approxEqual;
import std.parallelism;
import std.path;
import std.random: choice, dice, randomSample;
import std.range;
import std.stdio;
import std.uuid;

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
    double falsePositiveDetectionRate = 0.0;
    double grantApplicationCost = 0.0;
    string paramsList = "";
    bool syncFPRs;
    AwardPolicy policy = AwardPolicy.FPR;
    SelectionMethod selectionMethod = SelectionMethod.BEST_OF_TEN;

    auto helpInformation = getopt(
        args,
        std.getopt.config.passThrough,
        "nTrials", "Number of trials to run (default 10)", 
            &nTrials,
        "nIter", "Number of iterations",
            &N_ITER,
        "baseRate", "Base rate of true hypotheses (default 0.1)", 
            &baseRate,
        "awardAmount", "Amount given to grant-winning lab in a timestep (default 50)", 
            &awardAmount,
        "initialFalsePositiveRate", "False positive rate of all PIs at t=0 (default 0.05)", 
            &initialFalsePositiveRate,
        "fprMutationRate", "How often the false positive rate mutates (default 0.25)", 
            &fprMutationRate,
        "publishNegativeResultRate", "Rate that negative results are published (default 0.0)",
            &publishNegativeResultRate,
        "fprMutationMagnitude", "Std. dev. of the false positive mutations (default 0.01)",
            &fprMutationMagnitude,
        "falsePositiveDetectionRate", "Std. dev. of the false positive mutations (default 0.01)",
            &falsePositiveDetectionRate,
        "policy", "One of: RANDOM, PUBLICATIONS, FPR (default PUBLICATIONS)", 
            &policy,
        "selectionMethod", "One of: BEST_OF_TEN, WRIGHT_FISHER (default BEST_OF_TEN)", 
            &selectionMethod,
        "paramsList", "Comma-separated list of variable parameters, <POLICY>,<AWARD AMOUNT>,<PUB. NEG. RES. RATE>,<FALSE POS. DET. RATE>; e.g. \"FPR,5,0.5,0.9\"", 
            &paramsList,
        "syncFPRs", "Sync all agent FPR values at every synced timestep",
            &syncFPRs
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

    // comma-separated `paramsList` overrides options if provided.
    if (paramsList.length > 0) {

        auto params = paramsList.split(",");

        enforce(params.length == 4, 
                "paramsList must have four comma-separated values");

        policy = params[0].to!AwardPolicy;
        awardAmount = params[1].to!double;
        publishNegativeResultRate = params[2].to!double;
        falsePositiveDetectionRate = params[3].to!double;
    }

    defaultPoolThreads(totalCPUs);

    /* Set up simulation data storage */
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
        "publishNegativeResultRate": publishNegativeResultRate,
        "falsePositiveDetectionRate": falsePositiveDetectionRate
    ];

    data.meanFunds.length = nTrials;
    data.meanPublications.length = nTrials;
    data.sumFunds.length = nTrials;
    data.sumPublications.length = nTrials;
    data.medianFunds.length = nTrials;
    data.medianPublications.length = nTrials;

    data.falsePositiveRate.length = nTrials;
    data.falseDiscoveryRate.length = nTrials;

    data.agentFPRs.length = nTrials;

    /****** RUN TRIALS IN PARALLEL ******/
    foreach (trialIdx; parallel(nTrials.iota))
    {
        writefln("running %d of %d", trialIdx + 1, nTrials);

        TimeseriesData thisTrialData = simulation(
            policy, awardAmount, baseRate, initialFalsePositiveRate,
            fprMutationRate, fprMutationMagnitude, publishNegativeResultRate,
            falsePositiveDetectionRate, N_ITER, syncFPRs, selectionMethod
        );

        data.meanFunds[trialIdx] = thisTrialData.meanFunds;
        data.meanPublications[trialIdx] = thisTrialData.meanPublications;
        data.sumFunds[trialIdx] = thisTrialData.sumFunds;
        data.sumPublications[trialIdx] = thisTrialData.sumPublications;
        data.medianFunds[trialIdx] = thisTrialData.medianFunds;
        data.medianPublications[trialIdx] = thisTrialData.medianPublications;

        data.falsePositiveRate[trialIdx] = thisTrialData.falsePositiveRate;
        data.falseDiscoveryRate[trialIdx] = thisTrialData.falseDiscoveryRate;

        data.agentFPRs[trialIdx] = thisTrialData.agentFPRs;
    }

    writeDir
        .buildPath(randomUUID().to!string ~ ".json")
        .File("w")
        .write(
            data.serializeToJsonString
        );

    return 0;
}



const size_t N_PI = 100;
/* size_t N_ITER = 1e6.to!size_t; */
size_t SYNC_EVERY = 10000;
TimeseriesData simulation(AwardPolicy policy,
                double awardAmount, double baseRate, 
                double initialFalsePositiveRate, double fprMutationRate, 
                double fprMutationMagnitude, double publishNegativeResultRate,
                double falsePositiveDetectionRate, size_t N_ITER, 
                bool syncFPRs, SelectionMethod selectionMethod
                ) 
{
    PI[] pis; 

    // Set up random number generators.
    auto fprMutationAmountRange = 
        normalVar(0.0, fprMutationMagnitude).range;
    auto mutateFprNowRange = uniformVar(0.0, 1.0).range;

    /* t=0; initialize first generation of PIs */
    foreach (i; 0..N_PI)
    { 
        pis ~= new PI(baseRate);
        pis[i].funds = awardAmount;
        pis[i].falsePositiveRate = initialFalsePositiveRate;
        pis[i].publishNegativeResultRate = publishNegativeResultRate;
    }
    
    TimeseriesData data;
    data.meanFunds.length = N_ITER / SYNC_EVERY;
    data.meanPublications.length = N_ITER / SYNC_EVERY;
    data.sumFunds.length = N_ITER / SYNC_EVERY;
    data.sumPublications.length = N_ITER / SYNC_EVERY;
    data.medianFunds.length = N_ITER / SYNC_EVERY;
    data.medianPublications.length = N_ITER / SYNC_EVERY;

    data.falsePositiveRate.length = N_ITER / SYNC_EVERY;
    data.falseDiscoveryRate.length = N_ITER / SYNC_EVERY;

    data.agentFPRs.length = N_ITER / SYNC_EVERY;

    // Model iterations loop.
    foreach (iter; 0..N_ITER)
    {
        // PI's try to do research and publish it.
        pis.doScience(falsePositiveDetectionRate);
             
        // Select a PI to reproduce and a PI to die.
        /* pis.evolve(mutateFprNowRange, fprMutationAmountRange, */ 
        /*            awardAmount, fprMutationRate, SelectionMethod.BEST_OF_TEN); */
        pis.evolve(mutateFprNowRange, fprMutationAmountRange, 
                   awardAmount, fprMutationRate, selectionMethod);

        // After using the front element of random ranges, pop them off.
        mutateFprNowRange.popFront();
        fprMutationAmountRange.popFront();

        // Agency reviews "grant applications".
        pis.applyForGrants(awardAmount, policy);

        /* Sync model data for this timestep */
        if (iter % SYNC_EVERY == 0) 
        {
            size_t syncIdx = iter / SYNC_EVERY;

            /* Syncing funds information */
            data.sumFunds[syncIdx] = 
                pis.map!"a.funds".array.sum;
            data.meanFunds[syncIdx] = 
                pis.map!"a.funds".array.mean;

            // Calculate median funds.
            double[] funds = pis.map!"a.funds".array;
            funds.sort();
            data.medianFunds[syncIdx] = funds[$ / 2];

            /* Syncing publications information */
            data.sumPublications[syncIdx] = 
                pis.map!"a.publications".array.sum;
            data.meanPublications[syncIdx] = 
                pis.map!"a.publications".array.mean;

            double[] pubs = pis.map!"a.publications.to!double".array;
            pubs.sort();
            data.medianPublications[syncIdx] = pubs[$ / 2];

            /* Sync false positive rate */
            double[] agentFPRs = pis.map!"a.falsePositiveRate".array;
            data.falsePositiveRate[syncIdx] = agentFPRs.mean;
            if (syncFPRs)
                data.agentFPRs[syncIdx] = agentFPRs;

            /* Sync false discovery rate */
            data.falseDiscoveryRate[syncIdx] = pis.falseDiscoveryRate;
        }
        // Reset whether a PI published and if it was a false discovery.
        foreach (pi; pis)
        {
            pi.published = false;
            pi.falseDiscovery = false;
        }

    }

    return data;
}
unittest 
{
    AwardPolicy policy = AwardPolicy.FDR;
    double awardAmount = 50.0;
    double baseRate = 0.1;
    double initialFalsePositiveRate = 0.1;
    double fprMutationRate = 0.25;
    double fprMutationMagnitude = 0.01;
    double publishNegativeResultRate = 0.0;
    double falsePositiveDetectionRate = 0.0;
    size_t nIter = 10000;

    TimeseriesData simData = simulation(policy, awardAmount, baseRate, 
        initialFalsePositiveRate, fprMutationRate, fprMutationMagnitude,
        publishNegativeResultRate, falsePositiveDetectionRate, nIter,
        false, SelectionMethod.BEST_OF_TEN
    );
}


private double falseDiscoveryRate(PI[] pis)
{
    return 
        pis.map!(pi => pi.falseDiscovery.to!int).sum().to!double / 
            pis.map!(pi => pi.published.to!int).sum().to!double;
}
unittest
{
    PI[] pis;
    foreach (ii; 0..10)
        pis ~= new PI();
    
    foreach (ii; [0, 1, 2, 3])
        pis[ii].falseDiscovery = true;
    foreach (ii; [0, 1, 2, 3, 4, 5, 6, 7])
        pis[ii].published = true;

    assert(pis.falseDiscoveryRate == 0.5);
    
    foreach (ii; [2, 3])
        pis[ii].falseDiscovery = false;

    assert(pis.falseDiscoveryRate == 0.25);

    pis[0].falseDiscovery = false;
    pis[1].falseDiscovery = false;
    assert(pis.falseDiscoveryRate == 0.0);
}


private void doScience(PI[] pis, double falsePositiveDetectionRate) 
{
    foreach (pi; pis)
        pi.doScience(falsePositiveDetectionRate);
}


/****************** PRINCIPAL INVESTIGATOR *******************/
const static double SCIENCE_COST = 1.0;
const static double INIT_FUNDS = 10.0;
const double INIT_POWER = 0.8;
double BASE_RATE = 0.1;
const static double INIT_FALSE_POS_RATE = 0.05;
class PI {
    double funds = INIT_FUNDS;
    double power = INIT_POWER;
    double falsePositiveRate = INIT_FALSE_POS_RATE;
    double publishNegativeResultRate = 0.0;
    size_t publications = 0;
    size_t age = 0;
    float baseRate = 0.1;
    size_t falseDiscoveries = 0;
    bool falseDiscovery;
    bool published;
    UniformRange uniformRange;
    
    this() 
    {
        this.uniformRange = uniformVar(0.0, 1.0);
    }

    this(float baseRate)
    {
        this.uniformRange = uniformVar(0.0, 1.0);
        this.baseRate = baseRate;
    }

    this(UniformRange uniformRange) 
    {
        this.uniformRange = uniformRange;
    }

    this(UniformRange uniformRange, float baseRate) 
    {
        this.uniformRange = uniformRange;
        this.baseRate = baseRate;
    }

    public:
        void doScience(double falsePositiveDetectionRate=0.0)
        {
            // PI only does science if they have enough funds.
            if (this.funds >= SCIENCE_COST)
            {
                // Pay up.
                this.funds -= SCIENCE_COST;
                // The oracle says the hypothesis is True.
                if (hypothesisTrue())
                {
                    // PI found a true positive.
                    if (foundPositiveGivenTrue())
                    {
                        // True positive results always published.
                        this.publications += 1;
                        this.published = true;
                    }
                    // PI found a false negative.
                    else if (publishNegativeResult())
                    {
                        this.publications += 1;
                        this.published = true;
                        this.falseDiscovery = true;
                    }
                }
                else
                {
                    // PI found a false positive.
                    if (foundPositiveGivenFalse() && 
                        !falsePositiveDetected(falsePositiveDetectionRate))
                    {
                        this.publications += 1;
                        this.published = true;
                        this.falseDiscovery = true;
                    }
                    // PI found a true negative. 
                    else if (publishNegativeResult())
                    {
                        this.publications += 1;
                        this.published = true;
                    }
                }
            }
            ++age;
        }

        /**
         * Use when killing an old and birthing a new PI instead of deleting
         * and creating a new PI instance.
         */
        void reset()
        {
            this.funds = INIT_FUNDS;
            this.publications = 0;
            this.age = 0;
            this.falseDiscoveries = 0;
            this.falseDiscovery = false;
        }

    private:
        bool foundPositiveGivenTrue() 
        {
            bool ret = this.uniformRange.front < this.power;
            this.uniformRange.popFront();
            return ret;
        }

        bool foundPositiveGivenFalse()
        {
            bool ret = this.uniformRange.front < this.falsePositiveRate;
            this.uniformRange.popFront();
            return ret;
        }

        bool falsePositiveDetected(double falsePositiveDetectionRate)
        {

            bool ret = this.uniformRange.front < falsePositiveDetectionRate;
            this.uniformRange.popFront();
            return ret;
        }

        bool publishNegativeResult() 
        {
            bool ret = this.uniformRange.front < this.publishNegativeResultRate;
            this.uniformRange.popFront();
            return ret;
        }

        bool hypothesisTrue()
        {
            bool ret = this.uniformRange.front < this.baseRate;
            this.uniformRange.popFront();
            return ret;
        }
}
unittest {
    PI pi = new PI();

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

    // Now do science with peer review. Check that the number of publications
    // is within a range for different false positive 
    pi.reset();

    pi.power = 0.0;
    pi.publishNegativeResultRate = 0.0;
    pi.funds = 1000;

    // So every result is a false positive.
    BASE_RATE = 0.0; 
    pi.baseRate = 0.0; 
    pi.falsePositiveRate = 1.0;

    .writeln;
    foreach (fpDetectRate; [0.0, 0.2, 0.5, 0.8])
    {
        foreach (_; 1000.iota)
        {
            pi.doScience(fpDetectRate);
        }
        double expectedUpperLimit = 1000 - (1000*fpDetectRate) + (50 / (1 + fpDetectRate));
        double expectedLowerLimit = 1000 - (1000*fpDetectRate) - (50 / (1 + fpDetectRate));
        writefln(
            "FP Detect Rate: %.2f\nUpper: %.2f\nLower: %.2f\nPubs: %d\n",
            fpDetectRate, expectedUpperLimit, expectedLowerLimit, 
            pi.publications
        );
        if (fpDetectRate == 0.0)
        {
            assert(pi.publications == 1000);
        }
        else
        {
            assert(pi.publications < expectedUpperLimit);
            assert(pi.publications > expectedLowerLimit);
        }
        pi.publications = 0;
        pi.funds = 1000;
    }

    // Check that falseDiscovery can be cast to 1 if True, 0 if False
    assert([true, false, true, true].map!(a => a.to!size_t).reduce!"a+b" == 3);
}


/********* GRANT APPLICATIONS **********/
enum AwardPolicy {RANDOM, PUBLICATIONS, FPR, FDR}; 

enum SelectionMethod {WRIGHT_FISHER, BEST_OF_TEN};


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

        case AwardPolicy.FDR:

            applicants
                .minElement!"a.falseDiscoveries"
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
        double fprMutationRate,
        SelectionMethod selectionType
    )
{
    // "Kill" oldest PI.
    PI resetPI;
    resetPI = pis.randomSample(N_TO_DIE).maxElement!"a.age";
    auto sameAgePIs = pis.filter!(a => a.age == resetPI.age).array;
    if (sameAgePIs.length > 1)
    {
        resetPI = sameAgePIs.choice();    
    }
    resetPI.reset();
    resetPI.funds = initialFunds;
    
    // Reproduce the one with most funding from ten chosen at random.
    /* auto piToReproduce = pis.randomSample(N_TO_REPRODUCE) */
    /*                         .maxElement!"a.publications"; */
    
    PI piToReproduce;
    final switch (selectionType)
    {
        case SelectionMethod.BEST_OF_TEN:
            piToReproduce = pis.randomSample(N_TO_REPRODUCE)
                               .maxElement!"a.publications";
            break;

        case SelectionMethod.WRIGHT_FISHER:
            size_t reproducingIdx;
            /* foreach (i; 0..10) */
            /* { */
            /*     reproducingIdx = dice(pis.map!"a.publications"); */
            /*     writeln(reproducingIdx); */
            /* } */
            reproducingIdx = dice(pis.map!"a.publications");
            /* writeln("Pubs in evolve: ", pis.map!"a.publications"); */
            /* writeln("Repro idx: ", reproducingIdx); */
            piToReproduce = pis[reproducingIdx];
            break;
    }

    piToReproduce.reproduce(
        resetPI, mutateFprNowRange, 
        fprMutationAmountRange, fprMutationRate
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
        fprMutationRate, SelectionMethod.BEST_OF_TEN
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
            fprMutationRate, SelectionMethod.BEST_OF_TEN
        );
        mutateFprNowRange.popFront();
        fprMutationAmountRange.popFront();
    }
    assert(pis.map!"a.falsePositiveRate".sum != 1.0);

    size_t[] counts = 0UL.repeat().take(pis.length).array;
    foreach (int piIdx, ref PI pi; pis)
    {
        pi.publications = piIdx + 1;
        // Using FPR to identify which PI was selected; must be between 0 and 1.
        pi.falsePositiveRate = piIdx * 0.01;
    }
    foreach (trialIdx; 0..1e5)
    {
        double mu = 0.0;

        pis.evolve(
            mutateFprNowRange, fprMutationAmountRange, awardAmount,
            mu, SelectionMethod.WRIGHT_FISHER
        );

        // tell the selected one because it will have same falsePositiveRate
        // as another, but the other's publications won't match (they were
        // set equally above...this seems needlessly complicated but I'm 
        // going with it).
        size_t newIdx = pis.minIndex!"a.publications < b.publications";
        PI newPI = pis[newIdx];

        // New PI inherited selected PI's false positive rate, 
        // which is its index.
        size_t selectedIdx = (newPI.falsePositiveRate * 100).to!size_t;
        double selectedFPR = newPI.falsePositiveRate;
        counts[selectedIdx] += 1;

        pis[newIdx].falsePositiveRate = newIdx * 0.01;
        pis[newIdx].publications = newIdx + 1;
    }
    // Total publications (fitness) in population is 55.
    // Frequency each will be chosen to reproduce is piIdx+1 / 55. 
    auto expectedFrequencies = 1.iota(11).map!"a/55.0";
    auto calculatedFrequencies = counts.map!(a => a.to!float / 1e5);
    assert(approxEqual(expectedFrequencies, calculatedFrequencies, 1e-2, 1e-2),
           "\nexpected:\n" ~ expectedFrequencies.to!string ~ "\n\ncalculated:\n" 
           ~ calculatedFrequencies.to!string);
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


/* SYNC DATA STRUCTURES */
struct TrialsData
{
    Metadata metadata;
    double[][] meanFunds;
    double[][] meanPublications;
    double[][] sumFunds;
    double[][] sumPublications;
    double[][] medianFunds;
    double[][] medianPublications;
    double[][] falsePositiveRate;
    double[][] falseDiscoveryRate;
    double[][][] agentFPRs;
}


struct TimeseriesData
{
    double[] meanFunds;
    double[] meanPublications;
    double[] sumFunds;
    double[] sumPublications;
    double[] medianFunds;
    double[] medianPublications;
    double[] falsePositiveRate;
    double[] falseDiscoveryRate;
    double[][] agentFPRs;
}


struct Metadata
{
    double baseRate = 0.1;
    size_t syncEvery;
    string policy;
    double[string] parameters;
}
