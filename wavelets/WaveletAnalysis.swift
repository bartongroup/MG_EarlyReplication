//
//  WaveletAnalysis.swift
//  ReplicationWavelet
//
//  Created by Julian Blow on 09/11/2021.
//

import Foundation


func rickerWavelet(peakHeight: Double, numberOfDataPoints: Int) -> (xVals: [Double], yVals: [Double]) {
// Creates a Ricker Wavelet with x values between -1 and +1, and given peak height, with xDelta separation between consecutive x values. Returns an array of (x, y) values
    let maxXVal = 1.0                             // in principle this could be different; <1 truncates the arms, >1 is flat beyond 1
    let e = 2.71828
    var xVals = [Double]()
    var yVals = [Double]()
    guard numberOfDataPoints > 1 else {
        return (xVals, yVals)
    }
    let xDelta = 2 / Double(numberOfDataPoints - 1)
    let piSquared = Double.pi * Double.pi
//    for xVal in stride(from: -maxXVal, through: maxXVal, by: xDelta) {
    for dataPoint in 0..<numberOfDataPoints {
        let xVal = xDelta * Double(dataPoint) - maxXVal
        let squareProduct = piSquared * xVal * xVal
        let yVal = (1 - (2 * squareProduct)) * pow(e, -squareProduct) * peakHeight
        xVals.append(xVal)
        yVals.append(yVal)
    }
    return (xVals, yVals)
}

func waveletAnalysis(signalArray: [Double], wavelet : [Double]) -> [Double] {
// Returns an array of values where the wavelet array slides over the signalArray multiplying corresponding members and storing their sum.
// Values in the result array are centred on corresponding central values in the signalArray (ie the two arrays are aligned)
// At the edges of the signalArray, partial product sums are returned
    var result = [Double]()
    precondition(wavelet.count%2 == 1, "for proper alignment, wavelet array should have an odd number of values")
    let halfWaveletWidth = wavelet.count / 2                // truncates to get correct width
    for genomeIndex in 0..<(signalArray.count) {
        var productSum = 0.0
        for waveletIndex in 0..<wavelet.count {
            let genomeIndexOfProduct = genomeIndex - halfWaveletWidth + waveletIndex
            if (genomeIndexOfProduct >= 0) && (genomeIndexOfProduct < signalArray.count) {
                productSum += wavelet[waveletIndex] * signalArray[genomeIndexOfProduct]
            }
        }
        result.append(productSum)
    }
    return result
}


