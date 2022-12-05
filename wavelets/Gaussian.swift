//
//  Gaussian.swift
//  Replication_Wavelet
//
//  Created by Julian Blow on 28/03/2022.
//

import Foundation

struct Gaussian {
    let peakHeight: Double
    let centrePoint: Double
    let stDev: Double
     
    init(parameters: [Double]) {
        precondition(parameters.count == 3, "wrong number of parameters to initialise a Gaussian")
        peakHeight = parameters[0]
        centrePoint = parameters[1]
        stDev = parameters[2]
    }
    
    init(peakHeight: Double, centrePoint: Double, stDev: Double) {
        self.init(parameters: [peakHeight, centrePoint, stDev])
    }
    
    func valueFor(x: Double) -> Double {
        return peakHeight * exp(-(x - centrePoint) * (x - centrePoint) / (2 * stDev * stDev))
    }
    
    func inverseValueFor(y: Double) -> (Double, Double)? {
        // returns the two x values that would give the specified (y) output of the function
        // the function only has a defined output if the y value lies between zero and peak height
        if abs(y) > abs(peakHeight) { return nil }
        if y < 0 && peakHeight > 0 { return nil }
        if y > 0 && peakHeight < 0 { return nil }
        if peakHeight == 0 || y == 0 { return nil }
        let larger = centrePoint + sqrt(-2 * stDev * stDev * log(y / peakHeight))
        let smaller = centrePoint + centrePoint - larger
        return (smaller, larger)
    }
    
    func halfMaximumHeight() -> Double {
        return valueFor(x: centrePoint) / 2
    }
    
    func fullWidthAtPercentOfMaximumHeight(percent: Double) -> Double {                     // complicated if there is negative topWidth
        let xVals = inverseValueFor(y: percent * peakHeight / 100)
        if xVals == nil { return 0 }
        return xVals!.1 - xVals!.0
    }
    
    func fullWidthAtHalfMaximum() -> Double {              
        let peakHalfHeight = valueFor(x: centrePoint) / 2
        if peakHalfHeight == 0 { return 0 }
        let (smaller, larger) = inverseValueFor(y: peakHalfHeight)!       // guaranteed to give a defined result for this input
        return larger - smaller
    }
    
}
