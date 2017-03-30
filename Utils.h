#pragma once

// Restricts a value to be within the range [minimum, maximum].
template<typename T>
T clamp(T min, T max, T value)
{
    if (value <= min)
        return min;
    if (value >= max)
        return max;
    else
        return value;
}