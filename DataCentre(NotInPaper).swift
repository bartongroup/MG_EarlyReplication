//
//  DataCentre(NotInPaper).swift
//  Replication_Wavelet
//
//  Created by Julian Blow on 18/04/2022.
//

import Foundation

extension DataCentre {
    
    func optimalFTGaussianForSinglePeak(chromosome: String, startKb: Double, endKb: Double) -> [Double] {
        // queried with a single peak, returns the four parameters of the best-fitting flat-topped Gaussian found by Nelder-Mead
        // The sumSquaredError is the sum of the error (signal - prediction) squared
        // Nelder-Mead curve fitting returns if either sumSquaredError < (tolerance x signal sum) squared, or maxIterations have been performed
        // Nelder-Mead curve fitting requires decent guesses at the possible parameter ranges
        let replicationSignalArray = replicationDataSetDict[activeDataSetName]!.replicationSignalInRegion(chromosome: chromosome, startKb: startKb, endKb: endKb)
        let signalSum = replicationSignalArray.reduce(0, {$0 + $1} )
        let tolerance = 0.0001
        let toleratedError = (signalSum * tolerance) * (signalSum * tolerance)
        func fitFTGaussToReplicationSignal(parameters: [Double]) -> Double? {
              // takes the four parameters for flatToppedGaussian and evaluates the function against the data in replicationSignalArray
              // returns the r^2 error or nil if convergence has been achieved
              var totalError = 0.0
              let FTGauss = flatToppedGaussian(parameters: parameters)
              for index in 0..<replicationSignalArray.count {
                  let result = FTGauss.valueFor(x: Double(index))
                  totalError += (result - replicationSignalArray[index]) * (result - replicationSignalArray[index])
              }
              if totalError < toleratedError {
                  return nil                          // approximation good enough
              }
              else {
                  return totalError
              }
          }
                    // good initial parameter selection is critical, our best guesses and our estimated error in those guesses
                    // the parameters here only work when there is ~1 Mbp of signal surrounding the peak!!!
        let maxStepsPerIteration = 200                  // the number of Nelder-Mead steps, seems to be sufficient
        let maxIterations = 10                          // number of shuffles of the initial parameters tried - avoids the occasional misfit (~5% rate)
        let estPeakHeight = replicationSignalArray.max()!                   // guess low then go up from there
        let estPeakHeightError = estPeakHeight * 0.3                        // Nelder-Mead moves away from guess at these steps
        let estCentrePointBin = Double(replicationSignalArray.count) / 2
        let estCentrePointBinError = 20 / dataBinSizeKb                        // ± 20kb
        let estStDev = 50 / dataBinSizeKb
        let estStDevError = estStDev * 0.3
        let estTopWidth = 50 / dataBinSizeKb           // guess at around 50 kb and the method moves up from the guess
        let estTopWidthError = estTopWidth * 0.3
        let initialParams = [estPeakHeight, estCentrePointBin, estStDev, estTopWidth]
        let estimatedParameterError = [estPeakHeightError, estCentrePointBinError, estStDevError, estTopWidthError]
        let theNMFitter = NelderMead(modelFunction: fitFTGaussToReplicationSignal, initialParams: initialParams, estimatedParameterError: estimatedParameterError, maxStepsPerIteration: maxStepsPerIteration, maxIterations: maxIterations)
        theNMFitter.setParameterMinima([0,nil,0,0])
        let result = theNMFitter.NelderMeadFit()
        return result
    }
    
    
}
