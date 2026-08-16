"""
    CalculateStatistics(data::AbstractArray) -> Tuple

Calculate various statistical metrics for the given data.

# Arguments
- `data::AbstractArray`: An array of numerical data for which to calculate the statistics.

# Returns
- `Tuple`: A tuple containing the following statistical metrics:
  - `maximum_value::Number`: The maximum value in the data.
  - `index_max::CartesianIndex`: The index of the maximum value in the data.
  - `minimum_value::Number`: The minimum value in the data.
  - `index_min::CartesianIndex`: The index of the minimum value in the data.
  - `mean_value::Float64`: The mean (average) value of the data.
  - `median_value::Float64`: The median value of the data.
  - `sigma_value::Float64`: The standard deviation of the data.
  - `skewness_value::Float64`: The skewness of the data.
  - `kurtosis_value::Float64`: The kurtosis of the data.

# Description
This function calculates various statistical metrics for the provided data array. It computes the maximum and minimum values along with their indices, the mean, median, standard deviation (sigma), skewness, and kurtosis of the data. Non-finite values are ignored and extrema indices refer to the original array. An `ArgumentError` is thrown when no finite value remains.

# Example
```julia
# Example usage
data = randn(1000)
max_val, idx_max, min_val, idx_min, mean_val, median_val, sigma_val, skew_val, kurt_val = CalculateStatistics(data)
println("Maximum Value: ", max_val)
println("Index of Maximum Value: ", idx_max)
println("Minimum Value: ", min_val)
println("Index of Minimum Value: ", idx_min)
println("Mean Value: ", mean_val)
println("Median Value: ", median_val)
println("Standard Deviation: ", sigma_val)
println("Skewness: ", skew_val)
println("Kurtosis: ", kurt_val)
"""
function CalculateStatistics(data::AbstractArray)
    finite_indices = findall(isfinite, data)
    isempty(finite_indices) &&
        throw(ArgumentError("CalculateStatistics needs at least one finite value."))
    finite_data = data[finite_indices]

    maximum_value, local_index_max = findmax(finite_data)
    minimum_value, local_index_min = findmin(finite_data)
    index_max = finite_indices[local_index_max]
    index_min = finite_indices[local_index_min]

    mean_value = mean(finite_data)
    median_value = median(finite_data)
    sigma_value = std(finite_data)
    skewness_value = skewness(finite_data)
    kurtosis_value = kurtosis(finite_data)
    
    return maximum_value, index_max, minimum_value, index_min, mean_value, median_value, sigma_value, skewness_value, kurtosis_value
end

_stats_summary_values(stats) = Any[
    stats[1], Tuple(stats[2]), stats[3], Tuple(stats[4]),
    stats[5], stats[6], stats[7], stats[8], stats[9],
]

"""
    SummarizeStats(LOS::String, n::AbstractArray, T::AbstractArray, B1::AbstractArray, B2::AbstractArray, BLOS::AbstractArray, V1::AbstractArray, V2::AbstractArray, VLOS::AbstractArray) -> DataFrame

Calculate and summarize statistical metrics for given physical quantities along a specified line of sight (LOS).

# Arguments
- `LOS::String`: The line of sight direction ("x", "y", or "z").
- `n::AbstractArray`: Array of gas densities.
- `T::AbstractArray`: Array of gas temperatures.
- `B1::AbstractArray`: Array of the first magnetic field component.
- `B2::AbstractArray`: Array of the second magnetic field component.
- `BLOS::AbstractArray`: Array of the magnetic field component along the LOS.
- `V1::AbstractArray`: Array of the first velocity component.
- `V2::AbstractArray`: Array of the second velocity component.
- `VLOS::AbstractArray`: Array of the velocity component along the LOS.

# Returns
- `DataFrame`: A DataFrame containing maximum, extrema indices, mean, median,
  standard deviation, skewness, and kurtosis for each physical quantity.

# Description
This function calculates descriptive statistics for the physical quantities
along a specified line of sight. The results are organized into a DataFrame
with the appropriate component names for the selected LOS.

# Example
```julia
# Example usage
LOS = "z"
n = rand(100, 100, 100)
T = rand(100, 100, 100)
B1 = rand(100, 100, 100)
B2 = rand(100, 100, 100)
BLOS = rand(100, 100, 100)
V1 = rand(100, 100, 100)
V2 = rand(100, 100, 100)
VLOS = rand(100, 100, 100)

df = SummarizeStats(LOS, n, T, B1, B2, BLOS, V1, V2, VLOS)
println(df)
"""
function SummarizeStats(LOS, n, T, B1, B2, BLOS, V1, V2, VLOS)
    # Calculate statistics for each data series
    stats_n = CalculateStatistics(n)
    stats_T = CalculateStatistics(T)
    stats_B1 = CalculateStatistics(B1)
    stats_B2 = CalculateStatistics(B2)
    stats_BLOS = CalculateStatistics(BLOS)
    stats_V1 = CalculateStatistics(V1)
    stats_V2 = CalculateStatistics(V2)
    stats_VLOS = CalculateStatistics(VLOS)

    if LOS == "z"
        df = DataFrame(
            Quantity = ["Max", "IndMax", "Min", "IndMin", "Mean", "Median", "Std", "Skew", "Kurt"],
            n = _stats_summary_values(stats_n), T = _stats_summary_values(stats_T),
            Bx = _stats_summary_values(stats_B1), By = _stats_summary_values(stats_B2),
            Bz = _stats_summary_values(stats_BLOS), Vx = _stats_summary_values(stats_V1),
            Vy = _stats_summary_values(stats_V2), Vz = _stats_summary_values(stats_VLOS)
        )
    elseif LOS == "y"
        df = DataFrame(
            Quantity = ["Max", "IndMax", "Min", "IndMin", "Mean", "Median", "Std", "Skew", "Kurt"],
            n = _stats_summary_values(stats_n), T = _stats_summary_values(stats_T),
            Bx = _stats_summary_values(stats_B1), Bz = _stats_summary_values(stats_B2),
            By = _stats_summary_values(stats_BLOS), Vx = _stats_summary_values(stats_V1),
            Vz = _stats_summary_values(stats_V2), Vy = _stats_summary_values(stats_VLOS)
        )
    else
        df = DataFrame(
            Quantity = ["Max", "IndMax", "Min", "IndMin", "Mean", "Median", "Std", "Skew", "Kurt"],
            n = _stats_summary_values(stats_n), T = _stats_summary_values(stats_T),
            By = _stats_summary_values(stats_B1), Bz = _stats_summary_values(stats_B2),
            Bx = _stats_summary_values(stats_BLOS), Vy = _stats_summary_values(stats_V1),
            Vz = _stats_summary_values(stats_V2), Vx = _stats_summary_values(stats_VLOS)
        )
    end
    
    return df
end
