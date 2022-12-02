//
//  adjacentValueSimilarityMetric.swift
//  ReplicationWavelet
//

import Foundation

func adjacentValueSimilarityMetric(inputArray: [Double]) -> Double {
    /*
    Given an input array of n values, returns a metric which indicates the similarity between adjacent values in the input array.
    A result of 1 means that the adjacent values are ordered in a maximally similar way; 0 means random similarity, negative values indicate anti-similarity.
    These values are based on the n-1 absolute differences between adjacent values in the array (adjacencyDifferenceForArray).
    This is then compared with the absolute differences between adjacent values in a sorted version of the array (adjacencyDifferenceForSortedValues) and with the mean value of the absolute differences between adjacent values of all possible permutations of the values in the array (meanAdjacencyDifferenceForAllPermutations).
     The adjacencyDifferenceForSortedValues is equal to the absolute difference between the largest and smallest values in the array.
     The meanAdjacencyDifferenceForAllPermutations is equal to the average difference between all non-identical pairwise comparisons of the values in the input array.
     meanAdjacencyDifferenceForArray can therefore take values from adjacencyDifferenceForSortedValues (ordered) through meanAdjacencyDifferenceForAllPermutations
     (random) and higher (anti-ordered; to a complicated maximum value not explored further here).
     To project meanAdjacencyDifferenceForArray onto the output metric, we return:
     1 - ((adjacencyDifferenceForArray - adjacencyDifferenceForSortedValues) / (meanAdjacencyDifferenceForAllPermutations - adjacencyDifferenceForSortedValues))
     */
    if inputArray.count < 3 {
        return 0                           // no sense in arrays less than 3, so say they are random; 3 is moot as no difference between random and anti-similar
    }
    var adjacencyDifferenceSum = 0.0
    for index in 1..<inputArray.count {
        adjacencyDifferenceSum += abs(inputArray[index] - inputArray[index-1])
    }
    let adjacencyDifferenceForArray = adjacencyDifferenceSum / Double(inputArray.count - 1)
    let adjacencyDifferenceForSortedValues = (inputArray.max()! - inputArray.min()!) / Double(inputArray.count - 1)
    
    var sumOfAllPermDifferences = 0.0
    var numberOfAllPermDifferences = 0
    for index1 in 0..<inputArray.count {
        for index2 in (index1 + 1)..<inputArray.count {
            sumOfAllPermDifferences += abs(inputArray[index1] - inputArray[index2])
            numberOfAllPermDifferences += 1
        }
    }
    let meanAdjacencyDifferenceForAllPermutations = sumOfAllPermDifferences / Double(numberOfAllPermDifferences)
    let similarityMetric = 1 - ((adjacencyDifferenceForArray - adjacencyDifferenceForSortedValues) / (meanAdjacencyDifferenceForAllPermutations - adjacencyDifferenceForSortedValues))
    return similarityMetric
}
