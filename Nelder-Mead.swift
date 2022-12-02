//
//  Nelder-Mead.swift
//  Nelder-Mead
//
//  Created by Julian Blow on 10/12/2021.
//

import Foundation

class NelderMead {
    // Class to perform a NelderMead curve fit optimisation
    // modelFunction is the function called to test the parameters; its input is the parameter array, return is r^2 error or nil to terminate (fit achieved)
    // initialParams is an array of the initial guesses at the parameter values
    // estimatedParameterError is an estimate of how wide the search might be for each parameter and is used in setting up the initial simplex array
    // if there are minimum limits on the paramters, these can optionally be supplied
    // maxStepsPerIteration is the number of optimisation steps to be performed if the error is not reduced to ther required value
    // maxIterations is the maximum number of iterations; each iteration shuffles the initial paramters to avoid rare cases of sticking on a local minimum
    // the routine returns an array of the optimised parameter values
    struct parameterSet {
        let parameters: [Double]
        let error: Double
    }

    let modelFunction: ([Double]) -> Double?
    let initialParams: [Double]
    let numParams: Int
    let estimatedParameterError: [Double]
    var parameterMinima: [Double?]
    var maxStepsPerIteration:Int
    let maxIterations: Int
    
    init (modelFunction: @escaping ([Double]) -> Double?, initialParams: [Double], estimatedParameterError: [Double], maxStepsPerIteration:Int, maxIterations: Int) {
        self.modelFunction = modelFunction
        self.initialParams = initialParams
        numParams = initialParams.count
        self.estimatedParameterError = estimatedParameterError
        parameterMinima = Array(repeating: nil, count: numParams)
        self.maxStepsPerIteration = maxStepsPerIteration
        self.maxIterations = maxIterations
    }
    
    func setParameterMinima(_ minima: [Double?]) {
        parameterMinima = minima
    }
    
    func NelderMeadFit() -> [Double] {
        // performs the optimisation
        var bestParameterSet = [Double]()
        var bestError = Double.infinity
        for _ in 0..<maxIterations {
            var initialParameterArray = [[Double]]()
            for parameterIndex in 0..<numParams {            // make all the parameters for the n+1 parameter arrays
                var parameterValue = initialParams[parameterIndex]
                var setOfParameters = [Double]()
                for _ in 0...numParams {
                    setOfParameters.append(parameterValue)
                    parameterValue += estimatedParameterError[parameterIndex]
                }
                setOfParameters.shuffle()                    // so there is a random assortment
                initialParameterArray.append(setOfParameters)
            }
            
            var simplexArray = [parameterSet]()
            for simplexIndex in 0...numParams    {           // construct the n+1 simplexes using the shuffled set of parameters
                var parameterArray = [Double]()
                for parameterIndex in 0..<numParams {
                    parameterArray.append(initialParameterArray[parameterIndex][simplexIndex])
                }
                let result = modelFunction(parameterArray)
                if result == nil { return parameterArray }
                simplexArray.append(parameterSet(parameters: parameterArray, error: result!))
            }
            
            for _ in 0..<maxStepsPerIteration {
                simplexArray.sort(by: { $0.error < $1.error })
                                            // reflect worst through centroid of all but but worst
                var centroid = [Double]()                       // centroid of all but the worst parameter sets
                let simplexWithoutLast = simplexArray.dropLast()
                for parameterIndex in 0..<numParams {
                    centroid.append(simplexWithoutLast.reduce(0, { $0 + $1.parameters[parameterIndex] }) / Double(numParams))
                }
                var reflected = [Double]()               // reflection of the worst parameter set (simplexArray[numParams]) through the centroid
                var reflectionDelta = [Double]()           // vector for reflection of the worst parameter set (simplexArray[numParams]) through the centroid
                for parameterIndex in 0..<numParams {
                    var deltaComponent = centroid[parameterIndex] - simplexArray[numParams].parameters[parameterIndex]
                    if parameterMinima[parameterIndex] != nil {
                        if centroid[parameterIndex] + deltaComponent < parameterMinima[parameterIndex]! {
                            deltaComponent = parameterMinima[parameterIndex]! - centroid[parameterIndex]
                        }
                    }
                    reflectionDelta.append(deltaComponent)
                    reflected.append(centroid[parameterIndex] + deltaComponent)
                }
                let reflectedResult = modelFunction(reflected)
                if reflectedResult == nil { return reflected }
                if (reflectedResult! < simplexArray[numParams-1].error) && (reflectedResult! >= simplexArray[0].error) {   // if better than 2nd worst, but not  best
                    simplexArray[numParams] = parameterSet(parameters: reflected, error: reflectedResult!)                 // swap worst for reflected
                    continue
                }
                                            // possible extension of reflection
                if reflectedResult! < simplexArray[0].error {        // if reflected is better than best, extend the reflection (reflectionDelta) by gamma
                    var extended = [Double]()
                    for parameterIndex in 0..<numParams {
                        extended.append(reflected[parameterIndex] + reflectionDelta[parameterIndex])
                        if parameterMinima[parameterIndex] != nil {
                            if extended.last! < parameterMinima[parameterIndex]! {
                                extended = extended.dropLast()
                                extended.append(parameterMinima[parameterIndex]!)
                            }
                        }
                    }
                    let extendedResult = modelFunction(extended)
                    if extendedResult == nil { return extended }
                    if extendedResult! < reflectedResult! {
                        simplexArray[numParams] = parameterSet(parameters: extended, error: extendedResult!)                 // swap worst for extended
                        continue
                    } else {
                        simplexArray[numParams] = parameterSet(parameters: reflected, error: reflectedResult!)                  // swap worst for reflected
                        continue
                    }
                }
                                        // possible contraction, either inside or outside the original
                var contractedOutside = [Double]()
                for parameterIndex in 0..<numParams {
                    contractedOutside.append(centroid[parameterIndex] + (0.5 * reflectionDelta[parameterIndex]))
                    if parameterMinima[parameterIndex] != nil {
                        if contractedOutside.last! < parameterMinima[parameterIndex]! {
                            contractedOutside = contractedOutside.dropLast()
                            contractedOutside.append(parameterMinima[parameterIndex]!)
                        }
                    }
                }
                let contractedOutsideResult = modelFunction(contractedOutside)
                if contractedOutsideResult == nil { return contractedOutside }
                var contractedInside = [Double]()
                for parameterIndex in 0..<numParams {
                    contractedInside.append(centroid[parameterIndex] - (0.5 * reflectionDelta[parameterIndex]))
                }
                let contractedInsideResult = modelFunction(contractedInside)
                if contractedInsideResult == nil { return contractedInside }
                if (contractedOutsideResult! < contractedInsideResult!) && (contractedOutsideResult! < simplexArray[numParams].error) {
                    simplexArray[numParams] = parameterSet(parameters: contractedOutside, error: contractedOutsideResult!)         // swap worst for contractedOutside
                    continue
                }
                if  contractedInsideResult! < simplexArray[numParams].error {
                    simplexArray[numParams] = parameterSet(parameters: contractedInside, error: contractedInsideResult!)        // swap worst for contractedInside
                    continue
                }
                                                // final option is to shrink all points in the simplex towards the best
                for simplex in 1..<simplexArray.count {              // shrink all parameter sets except best at index 0
                    var parameterArray = [Double]()
                    for parameterIndex in 0..<numParams {
                        parameterArray.append((simplexArray[simplex].parameters[parameterIndex] + simplexArray[0].parameters[parameterIndex]) / 2)
                    }
                    let result = modelFunction(parameterArray)
                    if result == nil { return parameterArray }
                    simplexArray[simplex] = parameterSet(parameters: parameterArray, error: result!)
                }
            }
            simplexArray.sort(by: { $0.error < $1.error })
            if simplexArray[0].error < bestError {
                bestParameterSet = simplexArray[0].parameters
                bestError = simplexArray[0].error
            }
        }
        return bestParameterSet
    }
}
    
    



