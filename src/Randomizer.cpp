#include "Randomizer.h"

// Gets a random boolean value.
bool Randomizer::randomBool()
{
    return std::uniform_int_distribution<int>(0, 1)(rngEngine) == 1;
}

// Sets the seed for the generator. 
void Randomizer::setSeed(unsigned int newSeed)
{
    seed = newSeed;
    rngEngine.seed(newSeed);
}

// Gets the seed used by the generator to generate random numbers.
unsigned int Randomizer::getSeed()
{
    return seed;
}

unsigned int Randomizer::seed = std::random_device{}();

std::mt19937 Randomizer::rngEngine(seed);
