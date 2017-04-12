#pragma once

#include <random>
#include <type_traits>
#include <algorithm>

// Temporarily disable compiler warning C4544 for this file.
// ('unnamed-parameter': default template argument ignored on this template declaration)
#pragma warning(push)
#pragma warning(disable : 4544)

// A class for generating random numbers.
class Randomizer
{
public:
    // Gets a random integral value in the range [begin, end].
    template <typename T, std::enable_if_t<std::is_integral<T>::value, std::nullptr_t> = nullptr>
    static T random(T begin, T end);

    // Gets a random floating point value in the range [begin, end).
    template <typename T, std::enable_if_t<std::is_floating_point<T>::value, std::nullptr_t> = nullptr>
    static T random(T begin, T end);

    // Gets a random boolean value.
    static bool randomBool();

    // Gets a random element in the given container and returns a reference to it.
    template <typename C>
    static typename C::value_type& randomElement(C& container);

    // Gets a random element in the given container and returns a const reference to it.
    template <typename C>
    static const typename C::value_type& randomElement(const C& container);

    // Randomly reorders the elements in the given container.
    template <typename C>
    static C& shuffle(C& container);

    // Sets the seed for the generator. 
    static void setSeed(unsigned int seed);

    // Gets the seed used by the generator to generate random numbers.
    static unsigned int getSeed();

private:
    static unsigned int seed;

    static std::mt19937 rngEngine;
};

// Gets a random integral value in the range [begin, end].
template <typename T, std::enable_if_t<std::is_integral<T>::value, std::nullptr_t> = nullptr>
T Randomizer::random(T begin, T end)
{
    return std::uniform_int_distribution<T>(begin, end)(rngEngine);
}

// Gets a random floating point value in the range [begin, end).
template <typename T, std::enable_if_t<std::is_floating_point<T>::value, std::nullptr_t> = nullptr>
T Randomizer::random(T begin, T end)
{
    return std::uniform_real_distribution<T>(begin, end)(rngEngine);
}

// Gets a random element in the given container and returns a reference to it.
template <typename C>
typename C::value_type& Randomizer::randomElement(C& container)
{
    if (std::size(container) == 0)
        throw std::out_of_range("The given container is empty; it is not possible to select a random element.");

    auto it = std::begin(container);
    std::advance(it, random(0U, std::size(container) - 1U));
    return *it;
}

// Gets a random element in the given container and returns a const reference to it.
template <typename C>
const typename C::value_type& Randomizer::randomElement(const C& container)
{
    if (std::size(container) == 0)
        throw std::out_of_range("The given container is empty; it is not possible to select a random element.");

    auto it = std::begin(container);
    std::advance(it, random(0U, std::size(container) - 1U));
    return *it;
}

// Randomly reorders the elements in the given container.
template <typename C>
C& Randomizer::shuffle(C& container)
{
    std::shuffle(std::begin(container), std::end(container), rngEngine);
    return container;
}
#pragma warning(pop)
