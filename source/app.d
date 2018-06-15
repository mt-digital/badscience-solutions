import core.thread;

import std.algorithm;
import std.array;
import std.conv;
import std.getopt;
import std.parallelism;
import std.random: randomSample;
import std.range;
import std.stdio;

import mir.random;
import mir.random.engine;
import mir.random.variable;
import mir.math.common: fabs;
import mir.random.algorithm: range, RandomRange, sample;


alias RandomRange!(threadLocal!Random, UniformVariable!double) UniformRange;
alias RandomRange!(threadLocal!Random, NormalVariable!double) NormalRange;


void main () {
    writefln("Random enum: %s", AwardPolicy.RANDOM.to!string);
}



size_t N_PI = 100;
size_t N_ITER = 1_000_000;
size_t SYNC_EVERY = 2000;
void simulation(string writePath) 
{
    PI[] pis; 
    PI newPi;

    auto fprMutationAmountRange = normalVar(0.0, FALSE_POS_RATE_MUTATION_SD).range;
    auto mutateFprNowRange = uniformVar(0.0, 1.0).range;

    foreach (_; 0..N_PI)
        pis ~= new PI();

    foreach (iter; 0..N_ITER)
    {
        // PI's try to do research and publish it.
        pis.doScience;
             
        // Agency reviews "grant applications".
        pis.applyForGrants;

        // Select a PI to reproduce and a PI to die.
        pis.evolve(mutateFprNowRange, fprMutationAmountRange);

        if (iter % SYNC_EVERY == 0) {
            pis.writeState(writePath);
        }
    }
}


private void doScience(PI[] pis) 
{
    foreach (pi; pis)
        pi.doScience();
}


private void writeState(PI[] pis, string writeDir) 
{
    // TODO
}


/****************** PRINCIPAL INVESTIGATOR *******************/
const static double SCIENCE_COST = 1.0;
const static double INIT_FUNDS = 100.0;
const double INIT_POWER = 0.5;
const double BASE_RATE = 0.5;
const static double INIT_FALSE_POS_RATE = 0.6;
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
    writeln("Running tests!");
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
                .randomSample(1)
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
        NormalRange fprMutationAmountRange
    )
{
    // "Kill" oldest PI.
    auto resetPI = pis.randomSample(N_TO_DIE).maxElement!"a.age";
    resetPI.reset();
    
    // Reproduce the one with most funding from ten chosen at random.
    auto piToReproduce = pis.randomSample(N_TO_REPRODUCE)
                            .maxElement!"a.publications";

    piToReproduce.reproduce(
        resetPI, mutateFprNowRange, fprMutationAmountRange
    );
}


const double FALSE_POS_RATE_MUTATION_SD = 0.01;
const double FALSE_POS_RATE_MUTATION_PROBABILITY = 0.01;
private void reproduce(
        PI piToReproduce, 
        PI resetPI, 
        UniformRange mutateFprNowRange,
        NormalRange fprMutationAmountRange
    )
{
    double mutateAmt;
    bool mutate = mutateFprNowRange.front < FALSE_POS_RATE_MUTATION_PROBABILITY;
    mutateFprNowRange.popFront();
    if (mutate) 
    {
        // Get mutation amount and generate a new random FPR mutation value.
        mutateAmt = fprMutationAmountRange.front;
        fprMutationAmountRange.popFront();

        // Set FPR of the reset PI as the reproducing PI's FPR plus mutation.
        resetPI.falsePositiveRate = 
            piToReproduce.falsePositiveRate + mutateAmt;
    }
}
