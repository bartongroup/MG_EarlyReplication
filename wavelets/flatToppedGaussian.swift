//
//  flatToppedGaussian.swift
//  ReplicationWavelet
//
//  Created by Julian Blow on 11/12/2021.
//

import Foundation

struct flatToppedGaussian {
    // a gausssian curve with a flat section in the middle with width topWidth
    // if topwidth is < 0, omits the central topwidth portion of the gaussian
    let peakHeight: Double
    let centrePoint: Double
    let stDev: Double
    let topWidth: Double
     
    init(parameters: [Double]) {
        precondition(parameters.count == 4, "wrong number of parameters to initialise flatToppedGaussian")
        peakHeight = parameters[0]
        centrePoint = parameters[1]
        stDev = parameters[2]
        topWidth = parameters[3]
    }
    
    init(peakHeight: Double, centrePoint: Double, stDev: Double, topWidth: Double) {
        self.init(parameters: [peakHeight, centrePoint, stDev, topWidth])
    }
    
    func valueFor(x: Double) -> Double {
        let topLHS = centrePoint - (topWidth / 2)
        let topRHS =  topLHS + topWidth
        if x < topRHS && x > topLHS {
            return peakHeight
        }
        var mappedX: Double          // the x value to be used on the standard gaussion formula
        if x < centrePoint {
            mappedX = x + (topWidth / 2)
        }
        else {
            mappedX = x - (topWidth / 2)
        }
        return peakHeight * exp(-(mappedX - centrePoint) * (mappedX - centrePoint) / (2 * stDev * stDev))
    }
    
    func inverseValueFor(y: Double) -> (Double, Double)? {
        // returns the two x values that would give the specified (y) output of the function
        // the function only has a defined output if the y value lies between zero and peak height
        if abs(y) > abs(peakHeight) { return nil }
        if y < 0 && peakHeight > 0 { return nil }
        if y > 0 && peakHeight < 0 { return nil }
        if peakHeight == 0 || y == 0 { return nil }
        let larger = centrePoint + sqrt(-2 * stDev * stDev * log(y / peakHeight)) + (topWidth / 2)
        let smaller = centrePoint + centrePoint - larger
        return (smaller, larger)
    }
    
    func halfMaximumHeight() -> Double {
        if topWidth >= 0 {
            return peakHeight / 2
        }
        return valueFor(x: centrePoint) / 2
    }
    
    func fullWidthAtPercentOfMaximumHeight(percent: Double) -> Double {                     // complicated if there is negative topWidth
        let xVals = inverseValueFor(y: percent * peakHeight / 100)
        if xVals == nil { return 0 }
        return xVals!.1 - xVals!.0
    }
    
    func fullWidthAtHalfMaximum() -> Double {                     // complicated if there is negative topWidth
        if topWidth >= 0 {
            return abs(2 * sqrt(2 * log(2)) * stDev) + topWidth
        }
        let peakHalfHeight = valueFor(x: centrePoint) / 2
        if peakHalfHeight == 0 { return 0 }
        let (smaller, larger) = inverseValueFor(y: peakHalfHeight)!       // guaranteed to give a defined result for this input
        return larger - smaller
    }
    
}
