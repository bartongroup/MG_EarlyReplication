//
//  DataCentre.swift
//  ReplicationWavelet
//
//  Created by Julian Blow on 10/11/2021.
//

import Foundation

let miniumIsolatedPeakWidth = 200.0
let maximumIsolatedPeakWidth = 1200.0
let timepointDataSets = ["TM_10-40","TM_40-70", "TM_70-100", "TM_100-130"]

class DataCentre {
    let replicationDataSetDict: [String : EarlyRepDataSet]              // key = data set name
    let chromosomeSizeKb: [String: Double]                              // key = chromosome name
    let chromosomeNames: [String]
    let dataSetNames: [String]
    var activeDataSetName: String
    var activeDataSet: EarlyRepDataSet { return replicationDataSetDict[activeDataSetName]! }
    let dataBinSizeKb: Double
    let U2OSTimings: U2OSTimingDomains!
    
    init(path: String) {
        (replicationDataSetDict, dataSetNames, chromosomeSizeKb) = DataCentre.parseFile(path: path)
        chromosomeNames = Array(chromosomeSizeKb.keys)
        precondition(dataSetNames.count > 0, "failed to load any data sets")
        activeDataSetName = dataSetNames[0]
        let firstDataSet = replicationDataSetDict[activeDataSetName]!
        dataBinSizeKb = firstDataSet.replicationSignalsDict["chr1"]!.midPointsKb[1] - firstDataSet.replicationSignalsDict["chr1"]!.midPointsKb[0]
        U2OSTimings = U2OSTimingDomains()
    }
    
    func waveletWidthBinsFromKb(waveletPeakWidthKb: Double) -> Int {
    /* Returns the width of the wavelet as an even number of bins. The central 22.5% of Ricker wavelet is above zero, so we need to divide by 0.225, round to the nearest bin number and ensure there is an even number of bins for proper alignment */
        let desiredHalfWaveletWidthKb = waveletPeakWidthKb / (0.225 * 2)            // central 22.5% of Ricker wavelet is above zero
        return 2 * Int(round(desiredHalfWaveletWidthKb / dataBinSizeKb))            // doubling an Int ensures even number
    }
    
    func binMidpointKbp(binNumber: Int) -> Double {
        /* Returns the midpoint of the nth bin in kb */
        return (Double(binNumber) * dataBinSizeKb) + (dataBinSizeKb)
    }
    
    func analyseRegionWithWavelet(chromosome: String, startKb: Double, endKb: Double, waveletPeakWidthKb: Double) -> (genomePositionsKb: [Double], waveletResult: [Double])?  {
        // Analyses the specified region with a wavelet with the specified peak width
        // Returns a tuple with the midpoints of the bins in the region (in kb) and an array of the corresponding wavelet analysis results
        let waveletWidthBins = waveletWidthBinsFromKb(waveletPeakWidthKb: waveletPeakWidthKb)
        let waveletResultChrom = replicationDataSetDict[activeDataSetName]!.waveletAnalysisResults(forWaveletWidthBins: waveletWidthBins, chromosome: chromosome)
        let timingsChrom = replicationDataSetDict[activeDataSetName]!.replicationSignalsDict[chromosome]!
        
        let startBin = Int(startKb / dataBinSizeKb)
        guard startBin >= 0 else { return nil }
        guard startBin < waveletResultChrom.count - 1 else { return nil }
        let endBin = Int(endKb / dataBinSizeKb)
        guard endBin > startBin else { return nil }
        guard endBin < waveletResultChrom.count else { return nil }
        let waveletGenomePositions = Array(timingsChrom.midPointsKb[startBin...endBin])
        let waveletResultPortion = Array(waveletResultChrom[startBin...endBin])
        return (genomePositionsKb: waveletGenomePositions, waveletResult : waveletResultPortion)
    }
    
    func waveletPeaksInRegion(chromosome: String, startKb: Double, endKb: Double, waveletPeakWidthKb: Double) -> [waveletPeak]  {
        // Analyses the specified region with a wavelet with the specified peak width and returns details of each wavelet peaks in it
        let waveletWidthBins = waveletWidthBinsFromKb(waveletPeakWidthKb: waveletPeakWidthKb)
        let waveletPeaksArray = replicationDataSetDict[activeDataSetName]!.waveletPeaks(forWaveletWidthBins: waveletWidthBins)
        var peaksInRegion = [waveletPeak]()
        guard chromosomeSizeKb[chromosome] != nil else {
            return peaksInRegion
        }
        for peak in waveletPeaksArray {
            if (peak.chromosome == chromosome) && (peak.peakPositionKb >= startKb) && (peak.peakPositionKb <= endKb) {
                peaksInRegion.append(peak)
            }
        }
        return peaksInRegion
    }
    
    func waveletPeaksInTimingDomains(timingDomainArray: [TimingDomainDescription], waveletPeakWidthKb: Double) -> [waveletPeak] {
        // Returns an array of all wavelet peaks identified in the specified timing domains using a wavelet peak width as specified
        var peaksInRegion = [waveletPeak]()
        for timingDomain in timingDomainArray {
            peaksInRegion += waveletPeaksInRegion(chromosome: timingDomain.chromosome, startKb: (Double(timingDomain.start) / 1000), endKb: (Double(timingDomain.end) / 1000), waveletPeakWidthKb: waveletPeakWidthKb)
        }
        return peaksInRegion
    }
    
    func isolatedPeaksInTimingDomains(timingDomainArray: [TimingDomainDescription], waveletPeakWidthKb: Double, minimumHeight: Double) -> [waveletPeak] {
        // Returns an array of all wavelet peaks identified in the specified timing domains using a wavelet peak width as specified and whose wavelet peak height
        // is greater than minimumHeight and that are 'isolated'.
        // Isolated means that the total replication signal for a distance of wavelet peak width on either side of the peak is less than the specified
        // percentFlankingSignal of the signal under the wavelet peak
        let percentFlankingSignal = 25.0
        var peaksInTimingDomains = waveletPeaksInTimingDomains(timingDomainArray: timingDomainArray, waveletPeakWidthKb: waveletPeakWidthKb)
        peaksInTimingDomains = peaksInTimingDomains.filter { $0.peakHeight >= minimumHeight }
        var isolatedPeaks = [waveletPeak]()
        for peak in peaksInTimingDomains {
            let leftFlankStartKb = peak.startKb - waveletPeakWidthKb
            let rightFlankEndKb = peak.endKb + waveletPeakWidthKb
            let flankPlusPeakSignal = activeDataSet.replicationSignalInRegion(chromosome: peak.chromosome, startKb: leftFlankStartKb, endKb: rightFlankEndKb)
            let peakSignal = peak.totalReplicationSignal
            let flankingSignal = flankPlusPeakSignal.reduce(0, {$0 + $1}) - peakSignal
            if peakSignal * percentFlankingSignal > flankingSignal * 100 {
                isolatedPeaks.append(peak)
            }
        }
        return isolatedPeaks
    }
    
    func adjacentWaveletPeaks(waveletPeakWidthKb: Double, maximumDistanceBetweenPeaksKb: Double) -> [[waveletPeak]] {
        // Returns a list of  groups of peaks considered to be 'adjacent' - they are in the same timing domain and their peaks are less than the specified
        // maxmimum apart (maximumDistanceBetweenPeaksKb)
        let (earlyTDs, _) = U2OSTimings.allDomains()
        var groupsOfAdjacentPeaks = [[waveletPeak]]()
        for anEarlyTD in earlyTDs {
            var aGroupOfPeaks = [waveletPeak]()
            let peaksInTD = waveletPeaksInRegion(chromosome: anEarlyTD.chromosome, startKb: Double(anEarlyTD.start)/1000, endKb: Double(anEarlyTD.end)/1000, waveletPeakWidthKb: waveletPeakWidthKb)
            var previousPeakPositionKb = -(maximumDistanceBetweenPeaksKb * 2)        //  for new TD, set the previous peak so far away we must start a new group
            for aPeak in peaksInTD {
                if aPeak.peakPositionKb > previousPeakPositionKb + maximumDistanceBetweenPeaksKb {
                    groupsOfAdjacentPeaks.append(aGroupOfPeaks)
                    aGroupOfPeaks.removeAll()
                }
                aGroupOfPeaks.append(aPeak)
                previousPeakPositionKb = aPeak.peakPositionKb
            }
            groupsOfAdjacentPeaks.append(aGroupOfPeaks)                         // close down group at end of chromosome
        }
        groupsOfAdjacentPeaks = groupsOfAdjacentPeaks.filter({ $0.count > 0 })
        return groupsOfAdjacentPeaks
    }
    
    func neighbouringPeakHeights(groupOfAdjacentWavelets: [[waveletPeak]]) -> [(selfHeight: Double, neighbourHeight: Double)] {
        // Takes as input an array containing groups of neighbouring wavelet peaks as produced by adjacentWaveletPeaks() and for each peak creates one or two tuples
        // (depending on how many neighbours it has) with its own maximum replication signal and that of its neighbour
        var neighbourArray = [(selfHeight: Double, neighbourHeight: Double)]()
        for aGroupOfWavelets in groupOfAdjacentWavelets {
            if aGroupOfWavelets.count < 2 {
                continue
            }
            var arrayOfMaxRepSignals = [Double]()
            for peakIndex in 0..<aGroupOfWavelets.count {
                let thePeak = aGroupOfWavelets[peakIndex]
                let halfPeakWidth = Double(thePeak.waveletWidthBins) * dataBinSizeKb / 8       // the peak is a quarter of the wavelet, and we are taking half that
                let (_, maxRepSignal, _, _) = replicationDataSetDict[activeDataSetName]!.replicationMetricsInRegion(chromosome: thePeak.chromosome, startKb: thePeak.peakPositionKb-halfPeakWidth, endKb: thePeak.peakPositionKb+halfPeakWidth)
                arrayOfMaxRepSignals.append(maxRepSignal)
            }
            for peakIndex in 1..<arrayOfMaxRepSignals.count {
                neighbourArray.append((selfHeight: arrayOfMaxRepSignals[peakIndex-1], neighbourHeight: arrayOfMaxRepSignals[peakIndex]))
            }
        }
        return neighbourArray
    }
    
    func sizeSequenceOfAdjacentPeaks(waveletPeakWidthKb: Double, maximumDistanceBetweenPeaksKb: Double) -> [[Int]] {
        // Returns the size order of groups of peaks considered to be 'adjacent' - they are in the same timing domain and their peaks are less than the specified
        // maxmimum apart. The size is order is determined by the maximum replication signal in the peak.
        // The return array lists the peaks in order, with 1 for the largest, 2 for the next largest and so on.
        // Because the direction the group is encountered is arbitray we reverse the order of the array if the largest peak was found in the second half.
        let groupsOfAdjacentPeaks = adjacentWaveletPeaks(waveletPeakWidthKb: waveletPeakWidthKb, maximumDistanceBetweenPeaksKb: maximumDistanceBetweenPeaksKb)
        var sizeRankingsOfAdjacentPeaks = [[Int]]()
        for aGroupOfPeaks in groupsOfAdjacentPeaks {
            var peakInformation = [(repSignal: Double, indexInGroup: Int)]()
            for peakIndex in 0..<aGroupOfPeaks.count {
                let thePeak = aGroupOfPeaks[peakIndex]
                let halfPeakWidth = Double(thePeak.waveletWidthBins) * dataBinSizeKb / 8       // the peak is a quarter of the wavelet, and we are taking half that
                let (_, maxRepSignal, _, _) = replicationDataSetDict[activeDataSetName]!.replicationMetricsInRegion(chromosome: thePeak.chromosome, startKb: thePeak.peakPositionKb-halfPeakWidth, endKb: thePeak.peakPositionKb+halfPeakWidth)
                peakInformation.append((repSignal: maxRepSignal, indexInGroup: peakIndex))
            }
            let peaksSortedOnRepSignal = peakInformation.sorted(by: {$0.repSignal > $1.repSignal})
            var sizeRankings = Array(repeating: 0, count: aGroupOfPeaks.count)
            for index in 0..<aGroupOfPeaks.count {
                sizeRankings[peaksSortedOnRepSignal[index].indexInGroup] = index + 1
            }
            if peaksSortedOnRepSignal[0].indexInGroup >= aGroupOfPeaks.count / 2 {
                sizeRankings.reverse()
            }
            sizeRankingsOfAdjacentPeaks.append(sizeRankings)
        }
        return sizeRankingsOfAdjacentPeaks
    }
    
   func similarityInRepSignalsBetweenAdjacentPeaks(groupOfAdjacentWavelets:  [[waveletPeak]]) -> [Double] {
    // Takes as input an array containing groups of neighbouring wavelet peaks as produced by adjacentWaveletPeaks(). It then compares the similarity in max
    // replication signal between adjacent peaks using adjacentValueSimilarityMetric(), then returns an array of the similarities; in order to
    // adequately account for the separations, the value returned from adjacentValueSimilarityMetric() is repeated n-1 times for each n values in the group
        var adjacentPeakSimilarity = [Double]()
        for aGroupOfWavelets in groupOfAdjacentWavelets {
            if aGroupOfWavelets.count < 3 {                      // need at least 3 peaks to have a concept of neighbour versus non-neighbour
                continue
            }
            var arrayOfMaxRepSignals = [Double]()
            for peakIndex in 0..<aGroupOfWavelets.count {
                let thePeak = aGroupOfWavelets[peakIndex]
                let halfPeakWidth = Double(thePeak.waveletWidthBins) * dataBinSizeKb / 8       // the peak is a quarter of the wavelet, and we are taking half that
                let (_, maxRepSignal, _, _) = replicationDataSetDict[activeDataSetName]!.replicationMetricsInRegion(chromosome: thePeak.chromosome, startKb: thePeak.peakPositionKb-halfPeakWidth, endKb: thePeak.peakPositionKb+halfPeakWidth)
                arrayOfMaxRepSignals.append(maxRepSignal)
            }
            let peakSimilarity = adjacentValueSimilarityMetric(inputArray: arrayOfMaxRepSignals)
            adjacentPeakSimilarity += Array(repeating: peakSimilarity, count: arrayOfMaxRepSignals.count-1)
        }
        return adjacentPeakSimilarity
    }
    
    func peakSeparations(waveletPeakWidthKb: Double, minimumPeakHeight: Double, inTimingZones timingZones: [TimingDomainDescription]) -> [Double] {
        // Returns the separation between peaks identified in the specified timing zones
        let waveletWidthBins = waveletWidthBinsFromKb(waveletPeakWidthKb: waveletPeakWidthKb)
        let allWaveletPeaks = replicationDataSetDict[activeDataSetName]!.waveletPeaks(forWaveletWidthBins: waveletWidthBins)
        var peakSeparations = [Double]()
        
        for timingZone in timingZones {
            var waveletPeaksToAnalyse = allWaveletPeaks.filter { $0.chromosome == timingZone.chromosome }
            waveletPeaksToAnalyse = waveletPeaksToAnalyse.filter { $0.peakPositionKb >= Double(timingZone.start) / 1000 }
            waveletPeaksToAnalyse = waveletPeaksToAnalyse.filter { $0.peakPositionKb <= Double(timingZone.end) / 1000 }
            waveletPeaksToAnalyse = waveletPeaksToAnalyse.filter { $0.peakHeight >= minimumPeakHeight }
            waveletPeaksToAnalyse.sort { $0.peakPositionKb < $1.peakPositionKb }
            var previousWaveletPeak: waveletPeak? = nil
            for peak in waveletPeaksToAnalyse {
                if previousWaveletPeak != nil {
                    peakSeparations.append(peak.peakPositionKb - previousWaveletPeak!.peakPositionKb)
                }
                previousWaveletPeak = peak
            }
        }
        return peakSeparations
    }
    
    func waveletPeakHeightsInRegion(chromosome: String, startKb: Double, endKb: Double, waveletPeakWidthKb: Double) -> [Double]  {
        // Returns an array of the wavelet peak heagiths in the specified region.
        let peaksInRegion = waveletPeaksInRegion(chromosome: chromosome, startKb: startKb, endKb: endKb, waveletPeakWidthKb: waveletPeakWidthKb)
        let peakHeightsInRegion = peaksInRegion.map { $0.peakHeight }
        return peakHeightsInRegion
    }
    
    func waveletPeakHeightsInTimingDomains(timingDomainArray: [TimingDomainDescription], waveletPeakWidthKb: Double) -> [Double]  {
        // Returns an array of the wavelet peak heagiths in the specified array of Timing Domains.
        let peaksInRegion = waveletPeaksInTimingDomains(timingDomainArray: timingDomainArray, waveletPeakWidthKb: waveletPeakWidthKb)
        let peakHeightsInRegion = peaksInRegion.map { $0.peakHeight }
        return peakHeightsInRegion
    }
    
    func peaksPerTimingDomain(waveletPeakWidthKb:Double, minimumPeakHeight: Double, inTimingDomains timingDomains: [TimingDomainDescription]) -> [(numberOfPeaks: Int, sizeOfTimingDomainKb: Double)] {
        // Returns the number of wavelet peaks found in each timing domain in a tuple with the size of the domain
        var peaksPerTDArray = [(numberOfPeaks: Int, sizeOfTimingDomainKb: Double)]()
        for timingDomain in timingDomains {
            let peaks = waveletPeaksInRegion(chromosome: timingDomain.chromosome, startKb: (Double(timingDomain.start) / 1000), endKb: (Double(timingDomain.end) / 1000), waveletPeakWidthKb: waveletPeakWidthKb)
            let peaksAboveCutoff = peaks.filter { $0.peakHeight >= minimumPeakHeight }
            peaksPerTDArray.append((numberOfPeaks: peaksAboveCutoff.count, sizeOfTimingDomainKb: Double(timingDomain.end - timingDomain.start) / 1000))
        }
        return peaksPerTDArray
    }
    
    func totalRepSignalPerTimingDomain(timingDomains: [TimingDomainDescription]) -> [Double] {
        // Returns the sum of all the replication signals in the specified timing domains
        var repSignalPerTDArray = [Double]()
        for timingDomain in timingDomains {
            let (totalRepSignal, _, _, _) = replicationDataSetDict[activeDataSetName]!.replicationMetricsInRegion(chromosome: timingDomain.chromosome, startKb: Double(timingDomain.start)/1000, endKb: Double(timingDomain.end)/1000)
            repSignalPerTDArray.append(totalRepSignal)
        }
        return repSignalPerTDArray
    }
    
    func peakGrowthWaveletAnalysis(isolatedPeaks: [peakPositionDescriptor], searchWaveletPeakWidthKb: Double) -> [[Double]?] {
        // Uses the time point series in timepointDataSets to find the optimal wavelet widths for the peaks passed in as 'isolated peaks'
        // A nil return for a peak indicates that one or more of the widths in the time series was outwith the minimum or maximum width for the analysis
        let initialActiveDataSet = activeDataSetName
        var peakWidthsOverTime = [[Double]?]()                      // top level if is each peak; for each peak, width at each time point
        for aPeak in isolatedPeaks {
            var regionStartKb = aPeak.startKb - searchWaveletPeakWidthKb
            if regionStartKb < 0 { regionStartKb = 0 }
            var regionEndKb = aPeak.endKb + searchWaveletPeakWidthKb
            if regionEndKb > chromosomeSizeKb[aPeak.chromosome]! { regionEndKb = chromosomeSizeKb[aPeak.chromosome]! }
            var isolatedPeakWidthsAtThisPosition = [Double]()
            for aDataSet in timepointDataSets {
                activeDataSetName = aDataSet
                let (optimalWaveletPeakWidthKb, _) = optimalWaveletWidthForSinglePeak(chromosome: aPeak.chromosome, startKb: regionStartKb, endKb: regionEndKb)                         // search our original peak width on either side of the peak to get the optimal width
                if optimalWaveletPeakWidthKb > miniumIsolatedPeakWidth && optimalWaveletPeakWidthKb < maximumIsolatedPeakWidth {
                    isolatedPeakWidthsAtThisPosition.append(optimalWaveletPeakWidthKb)
                }
            }
            if isolatedPeakWidthsAtThisPosition.count == timepointDataSets.count {
                peakWidthsOverTime.append(isolatedPeakWidthsAtThisPosition)                     // only record peak if valid measurement at each time point (ie not minimum or maximum widths
            }
            else {
                peakWidthsOverTime.append(nil)
            }
        }
        activeDataSetName = initialActiveDataSet
        return peakWidthsOverTime
    }
    
    func optimalWaveletWidthForSinglePeak(chromosome: String, startKb: Double, endKb: Double) -> (waveletWidthKb: Double, waveletPositionKb: Double) {
        // Queried with a single peak, returns the width of the best-fitting wavelet peak width in kb
        var currentWaveletMax = -Double.greatestFiniteMagnitude
        var currentWaveletMaxPositionKb = 0.0
        var currentOptimalWavelet = 0.0
        for waveletWidthKb in stride(from: miniumIsolatedPeakWidth, through: maximumIsolatedPeakWidth, by: 25.0) {
            let analysisResults = analyseRegionWithWavelet(chromosome: chromosome, startKb: startKb, endKb: endKb, waveletPeakWidthKb: waveletWidthKb)
            guard let (positionsKb, waveletResults) = analysisResults else {
                return (0, 0)
            }
            let resultMaxEnumerated = waveletResults.enumerated().max(by: {$0.element < $1.element })
            let thisWaveletMax = resultMaxEnumerated!.element           // the best value for current wavelet peak width
            if thisWaveletMax > currentWaveletMax {
                currentWaveletMax = thisWaveletMax
                currentOptimalWavelet = waveletWidthKb
                currentWaveletMaxPositionKb = positionsKb[resultMaxEnumerated!.offset]
            }
        }
        return (waveletWidthKb : currentOptimalWavelet, waveletPositionKb: currentWaveletMaxPositionKb)
    }
    
    func peakGrowthGaussianAnalysis(isolatedPeaks: [peakPositionDescriptor], searchWaveletPeakWidthKb: Double) -> [[Double]?] {
        // Uses the time point series in timepointDataSets to find the optimal gaussian widths for the peaks passed in as 'isolated peaks'
        // A nil return for a peak indicates that one or more of the widths in the time series was outwith the minimum or maximum width for the analysis
        let initialActiveDataSet = activeDataSetName
        var peakWidthsOverTime = [[Double]?]()                      // top level if is each peak; for each peak, width at each time point
        for aPeak in isolatedPeaks {
            var regionStartKb = aPeak.startKb - searchWaveletPeakWidthKb
            if regionStartKb < 0 { regionStartKb = 0 }
            var regionEndKb = aPeak.endKb + searchWaveletPeakWidthKb
            if regionEndKb > chromosomeSizeKb[aPeak.chromosome]! { regionEndKb = chromosomeSizeKb[aPeak.chromosome]! }
            var isolatedPeakWidthsAtThisPosition = [Double]()
            for aDataSet in timepointDataSets {
                activeDataSetName = aDataSet
                let gaussianParameters = optimalGaussianForSinglePeak(chromosome: aPeak.chromosome, startKb: regionStartKb, endKb: regionEndKb)                         // search our original peak width on either side of the peak to get the optimal width
                let optimalGaussian = Gaussian(parameters: gaussianParameters)
                let optimalWaveletPeakWidthKb = optimalGaussian.fullWidthAtHalfMaximum() * dataBinSizeKb
                if optimalWaveletPeakWidthKb > miniumIsolatedPeakWidth && optimalWaveletPeakWidthKb < maximumIsolatedPeakWidth {
                    isolatedPeakWidthsAtThisPosition.append(optimalWaveletPeakWidthKb)
                }
            }
            if isolatedPeakWidthsAtThisPosition.count == timepointDataSets.count {
                peakWidthsOverTime.append(isolatedPeakWidthsAtThisPosition)                     // only record peak if valid measurement at each time point (ie not minimum or maximum widths
            }
            else {
                peakWidthsOverTime.append(nil)
            }
        }
        activeDataSetName = initialActiveDataSet
        return peakWidthsOverTime
    }
    
    func optimalGaussianForSinglePeak(chromosome: String, startKb: Double, endKb: Double) -> [Double] {
        // Queried with a single peak, returns the three parameters (height, centre point and standard deviation) of the best-fitting Gaussian found by Nelder-Mead
        // The sumSquaredError is the sum of the error (signal - prediction) squared
        // Nelder-Mead curve fitting returns if either sumSquaredError < (tolerance x signal sum) squared, or maxIterations have been performed
        // Nelder-Mead curve fitting requires decent guesses at the possible parameter ranges
        let replicationSignalArray = replicationDataSetDict[activeDataSetName]!.replicationSignalInRegion(chromosome: chromosome, startKb: startKb, endKb: endKb)
        let signalSum = replicationSignalArray.reduce(0, {$0 + $1} )
        let tolerance = 0.001
        let toleratedError = (signalSum * tolerance) * (signalSum * tolerance)
        func fitGaussToReplicationSignal(parameters: [Double]) -> Double? {
              // takes the three parameters for a Gaussian and evaluates the function against the data in replicationSignalArray
              // returns the r^2 error or nil if convergence has been achieved
              var totalError = 0.0
              let fittedGauss = Gaussian(parameters: parameters)
              for index in 0..<replicationSignalArray.count {
                  let result = fittedGauss.valueFor(x: Double(index))
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
        let maxIterations = 20                          // number of shuffles of the initial parameters tried - avoids the occasional misfit (~5% rate)
        let estPeakHeight = replicationSignalArray.max()!                   // guess low then go up from there
        let estPeakHeightError = estPeakHeight * 0.3                        // Nelder-Mead moves away from guess at these steps
        let estCentrePointBin = Double(replicationSignalArray.count) / 2
        let estCentrePointBinError = 20 / dataBinSizeKb                        // ± 20kb
        let estStDev = 50 / dataBinSizeKb
        let estStDevError = estStDev * 0.3
        let initialParams = [estPeakHeight, estCentrePointBin, estStDev]
        let estimatedParameterError = [estPeakHeightError, estCentrePointBinError, estStDevError]
        let theNMFitter = NelderMead(modelFunction: fitGaussToReplicationSignal, initialParams: initialParams, estimatedParameterError: estimatedParameterError, maxStepsPerIteration: maxStepsPerIteration, maxIterations: maxIterations)
        theNMFitter.setParameterMinima([0,nil,0,0])
        let result = theNMFitter.NelderMeadFit()
        return result
    }
    
    func allPeaksAtPosition(chromosome: String, peakPositionKb: Double, positionErrorKb: Double, inDataSets dataSetNames: [String], waveletPeakWidthKb: Double) -> [waveletPeak?] {
        // Returns an array of wavelet peaks with a peak at the specified position plus or minus the specified error in all the named data sets
        // If no peak is found in the data set, nil is entered into the array
        var waveletPeaks = [waveletPeak?]()
        let waveletWidthBins = waveletWidthBinsFromKb(waveletPeakWidthKb: waveletPeakWidthKb)
        for dataSetName in dataSetNames {
            let replicationDataSet = replicationDataSetDict[dataSetName]
            if replicationDataSet == nil { continue }
            let peak = replicationDataSet!.waveletPeakAtPosition(chromosome: chromosome, peakPositionKb: peakPositionKb, maxPositionErrorKb: positionErrorKb, waveletWidthBins: waveletWidthBins)
            waveletPeaks.append(peak)
        }
        return waveletPeaks
    }
    
    func allValleys(waveletPeakWidthKb: Double, edgeOffsetDueToForkElongationKb: Double) -> [(leftPeak: waveletPeak, rightPeak: waveletPeak)] {
        // Identifies 'valleys' in the wavelet signal: any place within an early timing domain flanked by two wavelet peaks, where the distance between the
        // two flanking peak edges is greater than an edge offset edgeOffsetDueToForkElongationKb which is the distance a fork might move over 2 hours
        // The edges of the valley are the edges of the flanking peaks minus edgeOffsetDueToForkElongationKb on either side
        var valleyList = [(leftPeak: waveletPeak, rightPeak: waveletPeak)]()
        let (earlyTDs, _) = U2OSTimings.allDomains()
        for timingDomain in earlyTDs {
            let peaksInTD = waveletPeaksInRegion(chromosome: timingDomain.chromosome, startKb: (Double(timingDomain.start) / 1000), endKb: (Double(timingDomain.end) / 1000), waveletPeakWidthKb: waveletPeakWidthKb)
            if peaksInTD.count < 2 { continue }                 // need at least 2 peaks in TD for a valley
            for peakIndex in 1..<peaksInTD.count {
                let leftPeak = peaksInTD[peakIndex - 1]
                let rightPeak = peaksInTD[peakIndex]
                let valleyLeftEdge = leftPeak.endKb + edgeOffsetDueToForkElongationKb
                let valleyRightEdge = rightPeak.startKb - edgeOffsetDueToForkElongationKb
                if valleyLeftEdge >= valleyRightEdge { continue }               // in case peaks are too close together for there to be a valley
                let aValley = (leftPeak: leftPeak, rightPeak: rightPeak)
                valleyList.append(aValley)
            }
        }
        return valleyList
    }
    
    func totalRepSignalInRegion(chromosome: String, startKb: Double, endKb: Double) -> Double {
        let (total, _, _, _) = replicationDataSetDict[activeDataSetName]!.replicationMetricsInRegion(chromosome: chromosome, startKb: startKb, endKb: endKb)
        return total
    }
    
    func maxRepSignalInRegion (chromosome: String, startKb: Double, endKb: Double) -> Double {
        let (_, maximum, _, _) = replicationDataSetDict[activeDataSetName]!.replicationMetricsInRegion(chromosome: chromosome, startKb: startKb, endKb: endKb)
        return maximum
    }
    
    func minRepSignalInRegion (chromosome: String, startKb: Double, endKb: Double) -> Double {
        let (_, _, minimum, _) = replicationDataSetDict[activeDataSetName]!.replicationMetricsInRegion(chromosome: chromosome, startKb: startKb, endKb: endKb)
        return minimum
    }
    func meanRepSignalInRegion (chromosome: String, startKb: Double, endKb: Double) -> Double {
        let (_, _, _, mean) = replicationDataSetDict[activeDataSetName]!.replicationMetricsInRegion(chromosome: chromosome, startKb: startKb, endKb: endKb)
        return mean
    }
    
    func repSignalMetricsInRegion (chromosome: String, startKb: Double, endKb: Double) -> (total: Double, maximum: Double, minimum: Double, mean: Double) {
        return replicationDataSetDict[activeDataSetName]!.replicationMetricsInRegion(chromosome: chromosome, startKb: startKb, endKb: endKb)
    }
    
    class func parseFile(path: String) -> (replicationDataSetDict: [String : EarlyRepDataSet], header: [String], chromosomeSizeKb: [String: Double]) {
        // In parsing the file, we need to separate out all the columns which are separate data sets. The columns are: chrom name, start posn, end posn, midpoint Mbp, data sets....
        // On the first pass through the data, we simply extract the columns. We then go through each data set in turn and put into into the right collections
        let networkFile = try? String(contentsOfFile: path)
        var fileLines = networkFile!.components(separatedBy: "\n")
        let headerComponents = fileLines[0].components(separatedBy: "\t")
        fileLines.removeFirst()
        let numberOfDataSets = headerComponents.count - 4       // columns = chrom name, start posn, end posn, midpoint Mbp, data sets....
        let firstDataLineComponents = fileLines[0].components(separatedBy: "\t")
        let dataBinSizeKb = Double(firstDataLineComponents[2])! / 1000
        var parsedReplicationDataSetDict = [String : EarlyRepDataSet]()                           // key = data set name
        var parsedChromosomeSizeKb = [String: Double]()                                            // key = chromosome name
        var arrayOfMidPointsKb = [Double]()
        var arrayOfRepSignals = Array(repeating: [Double](), count: numberOfDataSets)           // an array for each data set
        var currentChromosome = ""
        var chromosomeName = ""
        var updatedChromosomeSize = 0.0
        for line in fileLines {
            let components = line.components(separatedBy: "\t")
            if currentChromosome != components[0] {                                 // found a new chromosome
                if updatedChromosomeSize > 0 {                                      // we have data so not the first time we go through
                    precondition(parsedReplicationDataSetDict[chromosomeName] == nil, "trying to add a chromosome that already exists")
                    parsedChromosomeSizeKb[chromosomeName] = updatedChromosomeSize / 1000
                    for dataSetNumber in 0..<numberOfDataSets {
                        let aDataSetName = headerComponents[dataSetNumber + 4]
                        let theReplicationSignals = replicationSignals(chromosomeName: chromosomeName, dataSetName: aDataSetName, midPointsKb: arrayOfMidPointsKb, repSignal: arrayOfRepSignals[dataSetNumber])
                        if parsedReplicationDataSetDict[aDataSetName] == nil {                  // first pass we need to create the dictionary entry
                            let aReplicationDataSet = EarlyRepDataSet(dataSetName: aDataSetName, replicationSignalsDict: [chromosomeName : theReplicationSignals], dataBinSizeKb: dataBinSizeKb)
                            parsedReplicationDataSetDict[aDataSetName] = aReplicationDataSet
                        }
                        else {
                            parsedReplicationDataSetDict[aDataSetName]!.replicationSignalsDict[chromosomeName] = theReplicationSignals
                        }
                    }
                }
                chromosomeName = "chr" + components[0]
                currentChromosome = components[0]
                arrayOfMidPointsKb = [Double]()
                arrayOfRepSignals = Array(repeating: [Double](), count: numberOfDataSets)       // empty all intermediate arrays
                updatedChromosomeSize = 0
                precondition(Int(components[1]) == 0, "chromosome doesn't start at zero")
                if currentChromosome == "" || currentChromosome == "Y" {                        // omit Y chromsoome, cells are female!
                    break
                }
            }
            precondition(updatedChromosomeSize == Double(components[1])!, "gaps in chromosome \(chromosomeName) at \(updatedChromosomeSize)")
            updatedChromosomeSize = Double(components[2])!
            arrayOfMidPointsKb.append(Double(components[3])! * 1000)                            // midpoints are in Mbp
            for dataSetNumber in 0..<numberOfDataSets {
                arrayOfRepSignals[dataSetNumber].append(Double(components[dataSetNumber + 4])!)
            }
        }
        if arrayOfMidPointsKb.count > 0 {                                            // collect last chromosome (not elegant!!!!)
            precondition(parsedReplicationDataSetDict[chromosomeName] == nil, "trying to add a chromosome that already exists")
            parsedChromosomeSizeKb[chromosomeName] = updatedChromosomeSize / 1000
            for dataSetNumber in 0..<numberOfDataSets {
                let aDataSetName = headerComponents[dataSetNumber + 4]
                let theReplicationSignals = replicationSignals(chromosomeName: chromosomeName, dataSetName: aDataSetName, midPointsKb: arrayOfMidPointsKb, repSignal: arrayOfRepSignals[dataSetNumber])
                parsedReplicationDataSetDict[aDataSetName]!.replicationSignalsDict[chromosomeName] = theReplicationSignals
            }
        }
        let dataSetNames = Array(headerComponents.dropFirst(4))
        return (parsedReplicationDataSetDict, dataSetNames, parsedChromosomeSizeKb)
    }
    
}
