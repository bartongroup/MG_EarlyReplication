//
//  EarlyRepDataSet.swift
//  ReplicationWavelet
//
//  Created by Julian Blow on 14/11/2021.
//

import Foundation
import AppKit

struct replicationSignals {
    let chromosomeName: String
    let dataSetName: String
    let midPointsKb: [Double]
    let repSignal: [Double]
}

protocol peakPositionDescriptor {
    var chromosome: String { get }
    var startKb: Double { get }
    var endKb: Double { get }
}

struct peakPosition: peakPositionDescriptor {
    var chromosome: String
    var startKb: Double
    var endKb: Double
}

struct waveletPeak: peakPositionDescriptor {
    let waveletWidthBins: Int
    let dataBinSizeKb: Double
    let chromosome: String
    let peakPositionKb: Double
    let peakHeight: Double
    let startKb: Double
    let endKb: Double
    let totalReplicationSignal: Double
    let maxReplicationSignal: Double
    let minReplicationSignal: Double
    let meanBinReplicationSignal: Double
}

class EarlyRepDataSet {
    let replicationDataSetName: String
    var replicationSignalsDict: [String : replicationSignals]           // String key = chromosome name
    let dataBinSizeKb: Double
    private var waveletAnalysisResults = [Int : [String: [Double]]]()           // Int key = waveletWidthBins, String = chromosome name, array is wavelet array
    private var waveletPeaks = [Int : [waveletPeak]]()                          // Int key = waveletWidthBins
    
    init(dataSetName: String, replicationSignalsDict: [String : replicationSignals], dataBinSizeKb: Double) {
        replicationDataSetName = dataSetName
        self.replicationSignalsDict = replicationSignalsDict
        self.dataBinSizeKb = dataBinSizeKb
    }
    
    func replicationSignalInRegion(chromosome: String, startKb: Double, endKb: Double) -> [Double] {
        // returns an array containing the early replication signals in the specified region
        guard replicationSignalsDict[chromosome] != nil else {
            return []
        }
        let replicationData = replicationSignalsDict[chromosome]!
        var startBinIndex = Int(startKb / dataBinSizeKb)
        if startBinIndex < 0 {
            startBinIndex = 0
        }
        var endBinIndex = Int(ceil(endKb / dataBinSizeKb))
        if endBinIndex >= replicationData.repSignal.count {
            endBinIndex = replicationData.repSignal.count - 1
        }
        return Array(replicationData.repSignal[startBinIndex...endBinIndex])
    }
    
    func replicationMetricsInRegion(chromosome: String, startKb: Double, endKb: Double) -> (total: Double, maximum: Double, minimum: Double, mean: Double) {
        // returns the total, maximum, minimum and mean early replication signal in the specified region
        let signalArray = replicationSignalInRegion(chromosome: chromosome, startKb: startKb, endKb: endKb)
        if signalArray.count == 0 {
            return (total: 0, maximum: 0, minimum: 0, mean: 0)
        }
        let repSignalSum = signalArray.reduce(0, {$0 + $1 })
        let repSignalMax = signalArray.max()!
        let repSignalMin = signalArray.min()!
        let repSignalMean = repSignalSum / Double(signalArray.count)
        return (total: repSignalSum, maximum: repSignalMax, minimum: repSignalMin, mean: repSignalMean)
    }
    
    func waveletAnalysisResults(forWaveletWidthBins numBins: Int, chromosome: String) -> [Double] {
        // returns the wavelet analysis of the specified chromosome and wavelet width as defined by the number of data bins
        if waveletAnalysisResults[numBins] == nil {
            analyseDataSetWithWavelet(waveletWidthBins: numBins)
        }
        let results = waveletAnalysisResults[numBins]!
        return results[chromosome]!
    }
    
    func waveletPeaks(forWaveletWidthBins numBins: Int) -> [waveletPeak] {
        // returns all the wavelet peaks identified using a wavelet width as defined by the number of data bins
        if waveletPeaks[numBins] == nil {
            findWaveletPeaks(waveletWidthBins: numBins)
        }
        return waveletPeaks[numBins]!
    }
    
    func heatmap(chromosome: String, waveletWidthBins: [Int]) -> [[Double]] {
        // returns a heatmap of the wavelet analysis of the specified chromosome using wavelet widths (in data bin sizes) as defined by the waveletWidthBins parameter
        var heatmap = [[Double]]()
        for waveletWidth in waveletWidthBins {
            heatmap.append(waveletAnalysisResults(forWaveletWidthBins: waveletWidth, chromosome: chromosome))
        }
        return heatmap
    }
    
    func waveletPeakAtPosition(chromosome: String, peakPositionKb: Double, maxPositionErrorKb: Double, waveletWidthBins numBins: Int) -> waveletPeak? {
        // returns a wavelet peak which has a peak at the specified position if it exists, or nil if there is no peak there
        // maxPositionErrorKb specifies the degree of error permissible in matching the waveletPeak postion with the target peakPositionKb
        // waveletWidthBins defines the width of the wavelet
        let waveletPeakArray = waveletPeaks(forWaveletWidthBins: numBins)
        let minimumPosition = peakPositionKb - maxPositionErrorKb
        let maximumPosition = peakPositionKb + maxPositionErrorKb
        for aWaveletPeak in waveletPeakArray {
            if aWaveletPeak.chromosome != chromosome { continue }
            if aWaveletPeak.peakPositionKb < minimumPosition { continue }
            if aWaveletPeak.peakPositionKb > maximumPosition { continue }
            return aWaveletPeak
        }
        return nil
    }
    
    private func analyseDataSetWithWavelet(waveletWidthBins: Int) {
        // performs a wavelet analysis on the data set using a wavelet width as defined by the number of data bins
        let wavelet = rickerWavelet(peakHeight: 1/Double(waveletWidthBins), numberOfDataPoints: waveletWidthBins + 1)     // +1 to even value so wavelet array can be zero-centred
        var waveletResults = [String : [Double]]()
        for (chromName, signals) in replicationSignalsDict {
            waveletResults[chromName] = waveletAnalysis(signalArray: signals.repSignal, wavelet: wavelet.yVals)
        }
        waveletAnalysisResults[waveletWidthBins] = waveletResults
        findWaveletPeaks(waveletWidthBins: waveletWidthBins)
    }
    
    private func findWaveletPeaks(waveletWidthBins: Int) {
        // Generates the required data for the waveletPeak structure from the wavelet analysis numbers stored in waveletAnalysisResults
        if waveletPeaks[waveletWidthBins] != nil {
            return
        }
        if waveletAnalysisResults[waveletWidthBins] == nil {
            analyseDataSetWithWavelet(waveletWidthBins: waveletWidthBins)
        }
        var peaksArray = [waveletPeak]()
        let analysisDataSet = waveletAnalysisResults[waveletWidthBins]!
        for (chromosome, waveletResultArray) in analysisDataSet {
            var inPeak = false
            var currentPeakHeight = 0.0
            var currentPeakIndex = 0
            var startPeakIndex = 0
            for index in 0..<waveletResultArray.count {
                if waveletResultArray[index] > 0 {
                    if !inPeak {                                    // new peak
                        currentPeakHeight = 0
                        startPeakIndex = index
                        inPeak = true
                    }
                    if currentPeakHeight < waveletResultArray[index] {
                        currentPeakHeight = waveletResultArray[index]
                        currentPeakIndex = index
                    }
            
                } else {
                    if inPeak {                                     // end of peak
                        let peakPositionKb = Double(currentPeakIndex) * dataBinSizeKb
                        let peakStartKb = Double(startPeakIndex) * dataBinSizeKb
                        let peakEndKb = Double(index) * dataBinSizeKb
                        let (totalReplicationSignal, maxReplicationSignal, minReplicationSignal, meanReplicationSignal) = replicationMetricsInRegion(chromosome: chromosome, startKb: peakStartKb, endKb: peakEndKb)
                        peaksArray.append(waveletPeak(waveletWidthBins: waveletWidthBins, dataBinSizeKb: dataBinSizeKb, chromosome: chromosome, peakPositionKb: peakPositionKb, peakHeight: currentPeakHeight, startKb: peakStartKb, endKb: peakEndKb, totalReplicationSignal: totalReplicationSignal, maxReplicationSignal: maxReplicationSignal, minReplicationSignal: minReplicationSignal, meanBinReplicationSignal: meanReplicationSignal))
                        inPeak = false
                    }
                }
            }
        }
        waveletPeaks[waveletWidthBins] = peaksArray
    }
    
}
