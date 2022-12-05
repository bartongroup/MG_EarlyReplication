//
//  AppDelegate.swift
//  ReplicationWavelet
//
//  Created by Julian Blow on 09/11/2021.


import Cocoa
import AppKit

@main

class AppDelegate: NSObject, NSApplicationDelegate {

    @IBOutlet var window: NSWindow!
    @IBOutlet var dataFileTextField: NSTextField!
    @IBOutlet var waveletWidthTextField: NSTextField!
    @IBOutlet var dataSetPopupButton: NSPopUpButton!
    @IBOutlet var percentRemovalTextField: NSTextField!
    @IBOutlet var chromosomeNameTextField: NSTextField!
    @IBOutlet var startPointTextField: NSTextField!
    @IBOutlet var endPointTextField: NSTextField!
    @IBOutlet var visualiseButton: NSButton!
    @IBOutlet var showWaveletsButton: NSButton!
    @IBOutlet var heatmapButton: NSButton!
    @IBOutlet var timingDomainsButton: NSButton!
    @IBOutlet var autoPeakGrowthWaveletButton: NSButton!
    @IBOutlet var manualPeakGrowthWaveletButton: NSButton!
    @IBOutlet var autoPeakGrowthGaussianButton: NSButton!
    @IBOutlet var manualPeakGaussianButton: NSButton!
    @IBOutlet var earlyVLateButton: NSButton!
    @IBOutlet var runMetricsButton: NSButton!
    @IBOutlet var peakSeparationsButton: NSButton!
    @IBOutlet var peakWaveletSweepButton: NSButton!
    @IBOutlet var percentileSweepButton: NSButton!
    @IBOutlet var peaksPerTDButton: NSButton!
    @IBOutlet var valleyFillingButton: NSButton!
    @IBOutlet var peakActivationSeqButton: NSButton!
    @IBOutlet var peakNeighbourSimilarityButton: NSButton!
    @IBOutlet var fitWaveletSinglePeakButton: NSButton!
    @IBOutlet var fitFTGaussSinglePeakButton: NSButton!
    @IBOutlet var fitGaussianSinglePeakButton: NSButton!
    
    var theDataCentre: DataCentre!
    var dataFileName: String?
    var dataSetString: String {
        return "file " + dataFileName! + ", data set " + theDataCentre.activeDataSetName
    }
    
    @IBAction func loadReplicationData(sender: NSButton) {
        let openPanel = NSOpenPanel()
        openPanel.message = "file containing early replication data"
        openPanel.runModal()
        let panelResult = openPanel.url
        if panelResult == nil { return }
        let dataPathString = panelResult!.path
        theDataCentre = DataCentre(path: dataPathString)
        guard theDataCentre.replicationDataSetDict.count > 0 else {
            return
        }
        visualiseButton.isEnabled = true
        runMetricsButton.isEnabled = true
        dataSetPopupButton.removeAllItems()
        dataSetPopupButton.addItems(withTitles: theDataCentre.dataSetNames)
        theDataCentre.activeDataSetName = dataSetPopupButton.selectedItem!.title
        dataSetPopupButton.isEnabled = true
        peakSeparationsButton.isEnabled = true
        percentileSweepButton.isEnabled = true
        peakWaveletSweepButton.isEnabled = true
        heatmapButton.isEnabled = true
        timingDomainsButton.isEnabled = true
        peaksPerTDButton.isEnabled = true
        earlyVLateButton.isEnabled = true
        autoPeakGrowthWaveletButton.isEnabled = true
        manualPeakGrowthWaveletButton.isEnabled = true
        autoPeakGrowthGaussianButton.isEnabled = true
        manualPeakGaussianButton.isEnabled = true
        valleyFillingButton.isEnabled = true
        peakActivationSeqButton.isEnabled = true
        peakNeighbourSimilarityButton.isEnabled = true
        fitWaveletSinglePeakButton.isEnabled = true
        fitFTGaussSinglePeakButton.isEnabled = true
        fitGaussianSinglePeakButton.isEnabled = true
        dataFileName = panelResult!.lastPathComponent
        dataFileTextField.stringValue = dataFileName!
    }
    
    @IBAction func changeDataSet(sender: NSPopUpButton) {
        let dataSet = sender.selectedItem!
        theDataCentre.activeDataSetName = dataSet.title
    }
    
    @IBAction func visualiseData(sender: NSButton) {
        if !checkChromosomeParameters() {
            return
        }
        let chromosomeName = chromosomeNameTextField.stringValue
        let waveletPeakWidthKb = waveletWidthTextField.doubleValue
        let waveletWidthBins = theDataCentre.waveletWidthBinsFromKb(waveletPeakWidthKb: waveletPeakWidthKb)         // get wavelet width in units of data bins
        let startPointKb = startPointTextField.doubleValue * 1000
        let startBin = Int(startPointKb / theDataCentre.dataBinSizeKb)
        let endPointKb = endPointTextField.doubleValue * 1000
        let endBin = Int(endPointKb / theDataCentre.dataBinSizeKb)
    
        let (_, waveletResult) = theDataCentre.analyseRegionWithWavelet(chromosome: chromosomeName, startKb: startPointKb, endKb: endPointKb, waveletPeakWidthKb: waveletPeakWidthKb)!
        let scaleFactorForDisplay = Double(waveletWidthBins) * 0.2                // scale the wavelet result so it is a bit below the replication signal
        let waveletResultScaled = waveletResult.map { $0 * scaleFactorForDisplay }
        let waveletPeakHeightCutoff = minimumPeakHeightCutoffFromLateData(waveletPeakWidthKb: waveletPeakWidthKb, percentileCutoff: percentRemovalTextField.doubleValue) * scaleFactorForDisplay
        let cutOffLine: [(Double, Double)] = [(theDataCentre.binMidpointKbp(binNumber: startBin)/1000, waveletPeakHeightCutoff), (theDataCentre.binMidpointKbp(binNumber: endBin)/1000, waveletPeakHeightCutoff)]
        let replicationPortion = Array(theDataCentre.activeDataSet.replicationSignalsDict[chromosomeName]!.repSignal[startBin...endBin])
        
        let resultGraph = JBLineGraph()
        resultGraph.createEquallySpacedYDataSetNamed("replication", withYData: replicationPortion, firstX: startPointKb/1000, XDelta: theDataCentre.dataBinSizeKb/1000)
        resultGraph.setPlotType(.barGraph, dataSetName: "replication")
        resultGraph.setPlotFillColor(NSColor(red: 0.4, green: 0.8, blue: 1, alpha: 1), dataSetName: "replication")
        resultGraph.setLineWidth(0, dataSetName: "replication")
        resultGraph.setBarLineWidth(0, dataSetName: "replication")
        if showWaveletsButton.state == .on {
            resultGraph.createEquallySpacedYDataSetNamed("wavelet", withYData: waveletResultScaled, firstX: startPointKb/1000, XDelta: theDataCentre.dataBinSizeKb/1000)
            resultGraph.setPlotSymbolType(JBGraphSymbol.none, dataSetName: "wavelet")
            resultGraph.setLineColor(NSColor(red: 1, green: 0, blue: 0, alpha: 1), dataSetName: "wavelet")
            resultGraph.createXYDataSetNamed("cutoff", fromXYTupleArray: cutOffLine)
            resultGraph.setPlotSymbolType(JBGraphSymbol.none, dataSetName: "cutoff")
            resultGraph.setLineDashType(lineDashType.fineHalfDash, dataSetName: "cutoff")
            resultGraph.graphTitle = String(format: "wavelet peak width %.0f kb cutoff %0.f%%\n%@", waveletPeakWidthKb, percentRemovalTextField.doubleValue, dataSetString)
        }
        else {
            resultGraph.graphTitle = "early replication " + dataSetString
        }
        resultGraph.theXGraphAxis.userDefinedAxisVariables.axisMin = Decimal(startPointTextField.doubleValue)
        resultGraph.theXGraphAxis.userDefinedAxisVariables.axisMax = Decimal(endPointTextField.doubleValue)
        resultGraph.theYGraphAxis.userDefinedAxisVariables.axisMin = Decimal(0)
        resultGraph.setXAxisTitle("chromosome " + chromosomeName + " (Mbp)")
        resultGraph.setYAxisTitle("replication signal")
        resultGraph.createGraph(displayGraph: true)
    }
    
    func checkChromosomeParameters() -> Bool {
        guard theDataCentre.replicationDataSetDict.count > 0 else { return false }
        let chromosomeName = chromosomeNameTextField.stringValue
        guard theDataCentre.chromosomeSizeKb[chromosomeName] != nil else {
            chromosomeNameTextField.stringValue = "chr1"
            return false
        }
        let startPointKb = startPointTextField.doubleValue * 1000
        let startBin = Int(startPointKb / theDataCentre.dataBinSizeKb)
        guard  startBin >= 0 else {
            startPointTextField.doubleValue = 0
            return false
        }
        let endPointKb = endPointTextField.doubleValue * 1000
        guard endPointKb <= Double(theDataCentre.chromosomeSizeKb[chromosomeName]!) else {
            endPointTextField.doubleValue = Double(theDataCentre.chromosomeSizeKb[chromosomeName]! - 0.01) / 1000           // trim a little for rounding
            return false
        }
        let endBin = Int(endPointKb / theDataCentre.dataBinSizeKb)
        guard endBin > 3  else {
            startPointTextField.doubleValue = theDataCentre.dataBinSizeKb * Double(endBin - 1) / 1000
            return false
        }
        return true
    }
    
    @IBAction func heatmap(sender: NSButton) {
        if !checkChromosomeParameters() {
            return
        }
        let chromosomeName = chromosomeNameTextField.stringValue
        let startPointKb = startPointTextField.doubleValue * 1000
        let endPointKb = endPointTextField.doubleValue * 1000
        let peakWidthsForAssay: [Double] = [50,70.7,100,141,200,283,400,566,800,1131,1600,2263,3200,4525,6400,9051,12800,18102,25600,36204,51200]    // doubling plus intermediate values of previous  x sqrt(2)
        var densityGrid = [[Double]]()
        var maxWaveletPeakWidth = 0.0
        for peakWidth in peakWidthsForAssay {
            if peakWidth > (endPointKb - startPointKb) {                    // no point testing wavelet much wider than data range
                break
            }
            maxWaveletPeakWidth = peakWidth
            let (_, waveletResult) = theDataCentre.analyseRegionWithWavelet(chromosome: chromosomeName, startKb: startPointKb, endKb: endPointKb, waveletPeakWidthKb: peakWidth)!
            let filteredResults = waveletResult.map { $0 < 0 ? 0 : $0 }                     // only use positive values
            densityGrid.append(filteredResults)
        }
        let theHeatmapGraph = JBHeatmapGraph()
        theHeatmapGraph.createDensityDataSetNamed(nil, densityGrid: densityGrid, XValStart: startPointTextField.doubleValue, XValEnd: endPointTextField.doubleValue, YValStart: peakWidthsForAssay[0], YValEnd: peakWidthsForAssay.last!)
        theHeatmapGraph.setXAxisTitle("chromosome " + chromosomeName + " (Mbp)")
        theHeatmapGraph.setYAxisTitle("positive wavelet signal")
        theHeatmapGraph.graphTitle = String(format: "wavelet analysis (peak widths %.0f kb - %.0f kb)\n%@", peakWidthsForAssay[0], maxWaveletPeakWidth, dataSetString)
        theHeatmapGraph.createGraph(displayGraph: true)
    }
    
    @IBAction func showTimingDomains(sender: NSButton) {
        if !checkChromosomeParameters() {
            return
        }
        let chromosomeName = chromosomeNameTextField.stringValue
        let startPointKb = startPointTextField.doubleValue * 1000
        let endPointKb = endPointTextField.doubleValue * 1000
        let timingDomainData = theDataCentre.U2OSTimings.timingDomainLevelsInRegion(chromosome: chromosomeName, startKb: startPointKb, endKb: endPointKb)
        let earlyDomains = timingDomainData.map { $0 < 0 ? 0 : $0 }
        let lateDomains = timingDomainData.map { $0 > 0 ? 0 : $0 }
        let timingDomainsGraph = JBLineGraph()
        timingDomainsGraph.createEquallySpacedYDataSetNamed("early", withYData: earlyDomains, firstX: startPointKb/1000, XDelta: Double(theDataCentre.U2OSTimings.binSizeBp)/1000000)
        timingDomainsGraph.setPlotType(.barGraph, dataSetName: "early")
        timingDomainsGraph.setPlotFillColor(NSColor(red: 0.0, green: 0.75, blue: 0.0, alpha: 1), dataSetName: "early")
        timingDomainsGraph.setLineWidth(0, dataSetName: "early")
        timingDomainsGraph.setBarLineWidth(0, dataSetName: "early")
        timingDomainsGraph.createEquallySpacedYDataSetNamed("late", withYData: lateDomains, firstX: startPointKb/1000, XDelta: Double(theDataCentre.U2OSTimings.binSizeBp)/1000000)
        timingDomainsGraph.setPlotType(.barGraph, dataSetName: "late")
        timingDomainsGraph.setLineWidth(0, dataSetName: "late")
        timingDomainsGraph.setBarLineWidth(0, dataSetName: "late")
        timingDomainsGraph.setPlotFillColor(NSColor(red: 0.75, green: 0.0, blue: 0.0, alpha: 1), dataSetName: "late")
        timingDomainsGraph.setXAxisTitle("chromosome " + chromosomeName + " (Mbp)")
        timingDomainsGraph.setYAxisTitle("timing domain signal")
        timingDomainsGraph.graphTitle = "timing domains \(chromosomeName) \(startPointKb/1000) Mbp - \(endPointKb/1000) Mbp\n" + dataSetString
        timingDomainsGraph.createGraph(displayGraph: true)
    }
    
    func minimumPeakHeightCutoffFromLateData(waveletPeakWidthKb: Double, percentileCutoff: Double) -> Double {
        // finds all peak heights in late replicating timing domains, and sets waveletPeakHeightCutoff as a user-defined percentile of this
        let (_, lateTDs) = theDataCentre.U2OSTimings.allDomains()
        var waveletPeaksInLateDNA = theDataCentre.waveletPeaksInTimingDomains(timingDomainArray: lateTDs, waveletPeakWidthKb: waveletPeakWidthKb)
        waveletPeaksInLateDNA.sort { $0.peakHeight < $1.peakHeight}
        let sizeCutoffIndex = Int(Double(waveletPeaksInLateDNA.count - 1) * percentileCutoff / 100)
        return waveletPeaksInLateDNA[sizeCutoffIndex].peakHeight
    }
    
    @IBAction func autoPeakGrowthWaveletAnalysis(sender: NSButton) {
        let searchWaveletPeakWidthKb = waveletWidthTextField.doubleValue
        let isolatedPeaks = getAutoIsolatedPeaks(searchWaveletPeakWidthKb: searchWaveletPeakWidthKb)
        peakGrowthWaveletAnalysis(isolatedPeaks: isolatedPeaks, searchWaveletPeakWidthKb: searchWaveletPeakWidthKb, peakType: "auto")
    }
    
    func getAutoIsolatedPeaks(searchWaveletPeakWidthKb: Double) -> [waveletPeak] {
        let (earlyTDs, _) = theDataCentre.U2OSTimings.allDomains()
        let waveletPeakHeightCutoff = minimumPeakHeightCutoffFromLateData(waveletPeakWidthKb: searchWaveletPeakWidthKb, percentileCutoff: percentRemovalTextField.doubleValue)
        var isolatedPeaks = theDataCentre.isolatedPeaksInTimingDomains(timingDomainArray: earlyTDs, waveletPeakWidthKb: searchWaveletPeakWidthKb, minimumHeight: waveletPeakHeightCutoff)
        isolatedPeaks.sort { if $0.chromosome != $1.chromosome { return $0.chromosome < $1.chromosome } else { return $0.peakPositionKb < $1.peakPositionKb } }
        return isolatedPeaks
    }
    
    func peakGrowthWaveletAnalysis(isolatedPeaks: [peakPositionDescriptor], searchWaveletPeakWidthKb: Double, peakType: String) {
        let peakWidthsOverTime = theDataCentre.peakGrowthWaveletAnalysis(isolatedPeaks: isolatedPeaks, searchWaveletPeakWidthKb: searchWaveletPeakWidthKb)
        var outputString = "chromo\t start\tend\t10_40\t40-70\t70-100\t100-130\t" + dataFileName! + "\n"
        var TM10Widths = [Double]()
        var TM40Widths = [Double]()
        var TM70Widths = [Double]()
        var TM100Widths = [Double]()
        var TM40Rates = [Double]()
        var TM70Rates = [Double]()
        var TM100Rates = [Double]()
        var overallRates = [Double]()
        for aPeakSeriesIndex in 0..<isolatedPeaks.count {
            let isolatedPeak = isolatedPeaks[aPeakSeriesIndex]
            outputString += String(format: "%@\t%.2f\t%.2f\t", isolatedPeak.chromosome, isolatedPeak.startKb/1000, isolatedPeak.endKb/1000)
            if let aPeakSeries = peakWidthsOverTime[aPeakSeriesIndex] {
                for optimalWaveletPeakWidthKb in aPeakSeries {
                    outputString += String(format:"%.0f\t", optimalWaveletPeakWidthKb)
                }
                TM10Widths.append(aPeakSeries[0])
                TM40Widths.append(aPeakSeries[1])
                TM70Widths.append(aPeakSeries[2])
                TM100Widths.append(aPeakSeries[3])
                TM40Rates.append(aPeakSeries[1] - aPeakSeries[0])
                TM70Rates.append(aPeakSeries[2] - aPeakSeries[1])
                TM100Rates.append(aPeakSeries[3] - aPeakSeries[2])
                overallRates.append(aPeakSeries[3] - aPeakSeries[0])
            }
            outputString += "\n"
        }
        
        let savePanel = NSSavePanel()
        savePanel.nameFieldStringValue = peakType+" isolated wavelet peaks.tsv"
        savePanel.runModal()
        let panelResult = savePanel.url
        if panelResult == nil { return }
        let resultFileURL = panelResult?.absoluteURL
        try? outputString.write(toFile: resultFileURL!.path, atomically: true, encoding: String.Encoding.utf8)
        
        let TM10Mean = TM10Widths.reduce(0, { $0 + $1 }) / Double(TM10Widths.count)
        let TM40Mean = TM40Widths.reduce(0, { $0 + $1 }) / Double(TM40Widths.count)
        let TM70Mean = TM70Widths.reduce(0, { $0 + $1 }) / Double(TM70Widths.count)
        let TM100Mean = TM100Widths.reduce(0, { $0 + $1 }) / Double(TM100Widths.count)
        let TM10Variance = TM10Widths.reduce(0, { $0 + ($1-TM10Mean)*($1-TM10Mean) }) / Double(TM10Widths.count)
        let TM10Stdv = sqrt(TM10Variance)
        let TM40Variance = TM40Widths.reduce(0, { $0 + ($1-TM40Mean)*($1-TM40Mean) }) / Double(TM40Widths.count)
        let TM40Stdv = sqrt(TM40Variance)
        let TM70Variance = TM70Widths.reduce(0, { $0 + ($1-TM70Mean)*($1-TM70Mean) }) / Double(TM70Widths.count)
        let TM70Stdv = sqrt(TM70Variance)
        let TM100Variance = TM100Widths.reduce(0, { $0 + ($1-TM100Mean)*($1-TM100Mean) }) / Double(TM100Widths.count)
        let TM100Stdv = sqrt(TM100Variance)

        let TM10DataName = String(format: "10-40 widths, mean %.0f, stdev %.0f", TM10Mean, TM10Stdv)
        let TM40DataName = String(format: "40-70 widths, mean %.0f, stdev %.0f", TM40Mean, TM40Stdv)
        let TM70DataName = String(format: "70-100 widths, mean %.0f, stdev %.0f", TM70Mean, TM70Stdv)
        let TM100DataName = String(format: "100-130 widths, mean %.0f, stdev %.0f", TM100Mean, TM100Stdv)
        
        let TMAllGraph = JBLineGraph()
        TMAllGraph.createFrequencyDataSetNamed(TM10DataName, Data: TM10Widths, binSize: 100)
        TMAllGraph.setPlotType(.lineGraph, dataSetName: TM10DataName)
        TMAllGraph.createFrequencyDataSetNamed(TM40DataName, Data: TM40Widths, binSize: 100)
        TMAllGraph.setPlotType(.lineGraph, dataSetName: TM40DataName)
        TMAllGraph.createFrequencyDataSetNamed(TM70DataName, Data: TM70Widths, binSize: 100)
        TMAllGraph.setPlotType(.lineGraph, dataSetName: TM70DataName)
        TMAllGraph.createFrequencyDataSetNamed(TM100DataName, Data: TM100Widths, binSize: 100)
        TMAllGraph.setPlotType(.lineGraph, dataSetName: TM100DataName)
        TMAllGraph.setXAxisTitle("wavelet width (kb)")
        TMAllGraph.setYAxisTitle("frequency")
        TMAllGraph.graphTitle = String(format: "%d %@ isolated peak widths\n%@", TM10Widths.count, peakType, dataSetString)
        TMAllGraph.drawLineNames = true
        TMAllGraph.createGraph(displayGraph: true)
        
        let TM40RateMean = TM40Rates.reduce(0, { $0 + $1 }) / Double(TM40Rates.count)
        let TM70RateMean = TM70Rates.reduce(0, { $0 + $1 }) / Double(TM70Rates.count)
        let TM100RateMean = TM100Rates.reduce(0, { $0 + $1 }) / Double(TM100Rates.count)
        let overallRateMean = overallRates.reduce(0, { $0 + $1 }) / Double(overallRates.count)
        let TM40RatesStDev = sqrt(TM40Rates.reduce(0, {$0 + (($1-TM40RateMean) * ($1-TM40RateMean))} ) / Double(TM40Rates.count))
        let TM70RatesStDev = sqrt(TM70Rates.reduce(0, {$0 + (($1-TM70RateMean) * ($1-TM70RateMean))} ) / Double(TM70Rates.count))
        let TM100RatesStDev = sqrt(TM100Rates.reduce(0, {$0 + (($1-TM100RateMean) * ($1-TM100RateMean))} ) / Double(TM100Rates.count))
        let overallRatesStDev = sqrt(overallRates.reduce(0, {$0 + (($1-overallRateMean) * ($1-overallRateMean))} ) / Double(overallRates.count))
        let TM40RatesName = String(format: "10-40 rates, mean %.0f, stdev %.0f", TM40RateMean, TM40RatesStDev)
        let TM70RatesName = String(format: "40-70 rates, mean %.0f, stdev %.0f", TM70RateMean, TM70RatesStDev)
        let TM100RatesName = String(format: "70-100 rates, mean %.0f, stdev %.0f", TM100RateMean, TM100RatesStDev)
        let overallRatesName = String(format: "overall rates, mean %.0f, stdev %.0f", overallRateMean, overallRatesStDev)
        
        let allRatesGraph = JBLineGraph()
        allRatesGraph.createFrequencyDataSetNamed(TM40RatesName, Data: TM40Rates, binSize: 50)
        allRatesGraph.setPlotType(.lineGraph, dataSetName: TM40RatesName)
        allRatesGraph.createFrequencyDataSetNamed(TM70RatesName, Data: TM70Rates, binSize: 50)
        allRatesGraph.setPlotType(.lineGraph, dataSetName: TM70RatesName)
        allRatesGraph.createFrequencyDataSetNamed(TM100RatesName, Data: TM100Rates, binSize: 50)
        allRatesGraph.setPlotType(.lineGraph, dataSetName: TM100RatesName)
        allRatesGraph.createFrequencyDataSetNamed(overallRatesName, Data: overallRates, binSize: 50)
        allRatesGraph.setPlotType(.lineGraph, dataSetName: overallRatesName)
        allRatesGraph.setXAxisTitle("width increase (kb)")
        allRatesGraph.setYAxisTitle("frequency")
        allRatesGraph.graphTitle = String(format: "%d %@ width rate increases\n%@", TM10Widths.count, peakType, dataSetString)
        allRatesGraph.drawLineNames = true
        allRatesGraph.createGraph(displayGraph: true)
    }
    
    @IBAction func autoPeakGrowthGaussAnalysis(sender: NSButton) {
        let searchWaveletPeakWidthKb = waveletWidthTextField.doubleValue
        let isolatedPeaks = getAutoIsolatedPeaks(searchWaveletPeakWidthKb: searchWaveletPeakWidthKb)
        peakGrowthGaussAnalysis(isolatedPeaks: isolatedPeaks, waveletPeakWidthKb: searchWaveletPeakWidthKb, peakType: "auto")
    }
    
    func peakGrowthGaussAnalysis(isolatedPeaks: [peakPositionDescriptor], waveletPeakWidthKb: Double, peakType: String) {
        let peakWidthsOverTime = theDataCentre.peakGrowthGaussianAnalysis(isolatedPeaks: isolatedPeaks, searchWaveletPeakWidthKb: waveletPeakWidthKb)
        var outputString = "chromo\t start\tend\t10_40\t40-70\t70-100\t100-130\t" + dataFileName! + "\n"
        var TM10Widths = [Double]()
        var TM40Widths = [Double]()
        var TM70Widths = [Double]()
        var TM100Widths = [Double]()
        var TM40Rates = [Double]()
        var TM70Rates = [Double]()
        var TM100Rates = [Double]()
        var overallRates = [Double]()
        for aPeakSeriesIndex in 0..<isolatedPeaks.count {
            let isolatedPeak = isolatedPeaks[aPeakSeriesIndex]
            outputString += String(format: "%@\t%.2f\t%.2f\t", isolatedPeak.chromosome, isolatedPeak.startKb/1000, isolatedPeak.endKb/1000)
            if let aPeakSeries = peakWidthsOverTime[aPeakSeriesIndex] {
                for optimalWaveletPeakWidthKb in aPeakSeries {
                    outputString += String(format:"%.0f\t", optimalWaveletPeakWidthKb)
                }
                TM10Widths.append(aPeakSeries[0])
                TM40Widths.append(aPeakSeries[1])
                TM70Widths.append(aPeakSeries[2])
                TM100Widths.append(aPeakSeries[3])
                TM40Rates.append(aPeakSeries[1] - aPeakSeries[0])
                TM70Rates.append(aPeakSeries[2] - aPeakSeries[1])
                TM100Rates.append(aPeakSeries[3] - aPeakSeries[2])
                overallRates.append(aPeakSeries[3] - aPeakSeries[0])
            }
            outputString += "\n"
        }
        
        let savePanel = NSSavePanel()
        savePanel.nameFieldStringValue = peakType+" isolated Gaussian peaks.tsv"
        savePanel.runModal()
        let panelResult = savePanel.url
        if panelResult == nil { return }
        let resultFileURL = panelResult?.absoluteURL
        try? outputString.write(toFile: resultFileURL!.path, atomically: true, encoding: String.Encoding.utf8)
        
        let TM10Mean = TM10Widths.reduce(0, { $0 + $1 }) / Double(TM10Widths.count)
        let TM40Mean = TM40Widths.reduce(0, { $0 + $1 }) / Double(TM40Widths.count)
        let TM70Mean = TM70Widths.reduce(0, { $0 + $1 }) / Double(TM70Widths.count)
        let TM100Mean = TM100Widths.reduce(0, { $0 + $1 }) / Double(TM100Widths.count)
        let TM10Variance = TM10Widths.reduce(0, { $0 + ($1-TM10Mean)*($1-TM10Mean) }) / Double(TM10Widths.count)
        let TM10Stdv = sqrt(TM10Variance)
        let TM40Variance = TM40Widths.reduce(0, { $0 + ($1-TM40Mean)*($1-TM40Mean) }) / Double(TM40Widths.count)
        let TM40Stdv = sqrt(TM40Variance)
        let TM70Variance = TM70Widths.reduce(0, { $0 + ($1-TM70Mean)*($1-TM70Mean) }) / Double(TM70Widths.count)
        let TM70Stdv = sqrt(TM70Variance)
        let TM100Variance = TM100Widths.reduce(0, { $0 + ($1-TM100Mean)*($1-TM100Mean) }) / Double(TM100Widths.count)
        let TM100Stdv = sqrt(TM100Variance)

        let TM10DataName = String(format: "10-40 widths, mean %.0f, stdev %.0f", TM10Mean, TM10Stdv)
        let TM40DataName = String(format: "40-70 widths, mean %.0f, stdev %.0f", TM40Mean, TM40Stdv)
        let TM70DataName = String(format: "70-100 widths, mean %.0f, stdev %.0f", TM70Mean, TM70Stdv)
        let TM100DataName = String(format: "100-130 widths, mean %.0f, stdev %.0f", TM100Mean, TM100Stdv)
        
        let TMAllGraph = JBLineGraph()
        TMAllGraph.createFrequencyDataSetNamed(TM10DataName, Data: TM10Widths, binSize: 100)
        TMAllGraph.setPlotType(.lineGraph, dataSetName: TM10DataName)
        TMAllGraph.createFrequencyDataSetNamed(TM40DataName, Data: TM40Widths, binSize: 100)
        TMAllGraph.setPlotType(.lineGraph, dataSetName: TM40DataName)
        TMAllGraph.createFrequencyDataSetNamed(TM70DataName, Data: TM70Widths, binSize: 100)
        TMAllGraph.setPlotType(.lineGraph, dataSetName: TM70DataName)
        TMAllGraph.createFrequencyDataSetNamed(TM100DataName, Data: TM100Widths, binSize: 100)
        TMAllGraph.setPlotType(.lineGraph, dataSetName: TM100DataName)
        TMAllGraph.setXAxisTitle("gaussian FWHM width (kb)")
        TMAllGraph.setYAxisTitle("frequency")
        TMAllGraph.graphTitle = String(format: "%d %@ isolated peak widths\n%@", TM10Widths.count, peakType, dataSetString)
        TMAllGraph.drawLineNames = true
        TMAllGraph.createGraph(displayGraph: true)
        
        let TM40RateMean = TM40Rates.reduce(0, { $0 + $1 }) / Double(TM40Rates.count)
        let TM70RateMean = TM70Rates.reduce(0, { $0 + $1 }) / Double(TM70Rates.count)
        let TM100RateMean = TM100Rates.reduce(0, { $0 + $1 }) / Double(TM100Rates.count)
        let overallRateMean = overallRates.reduce(0, { $0 + $1 }) / Double(overallRates.count)
        let TM40RatesStDev = sqrt(TM40Rates.reduce(0, {$0 + (($1-TM40RateMean) * ($1-TM40RateMean))} ) / Double(TM40Rates.count))
        let TM70RatesStDev = sqrt(TM70Rates.reduce(0, {$0 + (($1-TM70RateMean) * ($1-TM70RateMean))} ) / Double(TM70Rates.count))
        let TM100RatesStDev = sqrt(TM100Rates.reduce(0, {$0 + (($1-TM100RateMean) * ($1-TM100RateMean))} ) / Double(TM100Rates.count))
        let overallRatesStDev = sqrt(overallRates.reduce(0, {$0 + (($1-overallRateMean) * ($1-overallRateMean))} ) / Double(overallRates.count))
        let TM40RatesName = String(format: "10-40 rates, mean %.0f, stdev %.0f", TM40RateMean, TM40RatesStDev)
        let TM70RatesName = String(format: "40-70 rates, mean %.0f, stdev %.0f", TM70RateMean, TM70RatesStDev)
        let TM100RatesName = String(format: "70-100 rates, mean %.0f, stdev %.0f", TM100RateMean, TM100RatesStDev)
        let overallRatesName = String(format: "overall rates, mean %.0f, stdev %.0f", overallRateMean, overallRatesStDev)
        
        let allRatesGraph = JBLineGraph()
        allRatesGraph.createFrequencyDataSetNamed(TM40RatesName, Data: TM40Rates, binSize: 50)
        allRatesGraph.setPlotType(.lineGraph, dataSetName: TM40RatesName)
        allRatesGraph.createFrequencyDataSetNamed(TM70RatesName, Data: TM70Rates, binSize: 50)
        allRatesGraph.setPlotType(.lineGraph, dataSetName: TM70RatesName)
        allRatesGraph.createFrequencyDataSetNamed(TM100RatesName, Data: TM100Rates, binSize: 50)
        allRatesGraph.setPlotType(.lineGraph, dataSetName: TM100RatesName)
        allRatesGraph.createFrequencyDataSetNamed(overallRatesName, Data: overallRates, binSize: 50)
        allRatesGraph.setPlotType(.lineGraph, dataSetName: overallRatesName)
        allRatesGraph.setXAxisTitle("Gaussian width increase (kb)")
        allRatesGraph.setYAxisTitle("frequency")
        allRatesGraph.graphTitle = String(format: "%d %@ width rate increases\n%@", TM10Widths.count, peakType, dataSetString)
        allRatesGraph.drawLineNames = true
        allRatesGraph.createGraph(displayGraph: true)
    }
    
    @IBAction func manualPeakGrowthWaveletAnalysis(sender: NSButton) {
        let isolatedPeaks = getManualPeaks()
        let waveletPeakWidthKb = waveletWidthTextField.doubleValue
        peakGrowthWaveletAnalysis(isolatedPeaks: isolatedPeaks, searchWaveletPeakWidthKb: waveletPeakWidthKb, peakType: "manual")
    }
    
    func getManualPeaks() -> [peakPositionDescriptor] {
        // Loads chromosome name and the start point and end point (both in Mbp) of a single peak into a tuple */
//        var manualPeaks = [waveletPeak]()
        var manualPeaks = [peakPositionDescriptor]()
        let openPanel = NSOpenPanel()
        openPanel.message = "list of isloated peaks to analyse"
        openPanel.runModal()
        let panelResult = openPanel.url
        if panelResult == nil { return manualPeaks }
        let dataPathString = panelResult!.path
        let peaksFile = try? String(contentsOfFile: dataPathString)
        
  //      theDataCentre.activeDataSetName = "TM_10-40"
  //      let waveletPeakWidthKb = waveletWidthTextField.doubleValue
        let fileLines = peaksFile!.components(separatedBy: CharacterSet(charactersIn: "\n\r"))
        for aLine in fileLines {
            if aLine.hasPrefix("//") || aLine.isEmpty {                  // ignore comment lines starting //
                continue
            }
            let lineComponents = aLine.components(separatedBy: "\t")
            let aPeak = peakPosition(chromosome: lineComponents[0], startKb: Double(lineComponents[1])! * 1000, endKb: Double(lineComponents[2])! * 1000)
            manualPeaks.append(aPeak)
/*            let peaks = theDataCentre.waveletPeaksInRegion(chromosome: lineComponents[0], startKb: Double(lineComponents[1])! * 1000, endKb: Double(lineComponents[2])! * 1000, waveletPeakWidthKb: waveletPeakWidthKb)
            if peaks.count == 1 {
                manualPeaks.append(peaks[0])
            }
            else {
                print("\(lineComponents[0]) \(lineComponents[1]) - \(lineComponents[2])")
            } */
        }
        manualPeaks.sort { if $0.chromosome != $1.chromosome { return $0.chromosome < $1.chromosome } else { return $0.startKb < $1.startKb } }
//        theDataCentre.activeDataSetName = dataSetPopupButton.selectedItem!.title
        return manualPeaks
    }
    
    @IBAction func manualPeakGrowthGaussAnalysis(sender: NSButton) {
        let isolatedPeaks = getManualPeaks()
        let waveletPeakWidthKb = waveletWidthTextField.doubleValue
        peakGrowthGaussAnalysis(isolatedPeaks: isolatedPeaks, waveletPeakWidthKb: waveletPeakWidthKb, peakType: "manual")
    }
    
    @IBAction func singlePeakWaveletFit(sender: NSButton) {
        if !checkChromosomeParameters() {
            return
        }
        let chromosomeName = chromosomeNameTextField.stringValue
        let startPointKb = startPointTextField.doubleValue * 1000
        let startBin = Int(startPointKb / theDataCentre.dataBinSizeKb)
        let endPointKb = endPointTextField.doubleValue * 1000
        let endBin = Int(endPointKb / theDataCentre.dataBinSizeKb)
        let dataSetsToTest = ["TM_10-40", "TM_40-70", "TM_70-100", "TM_100-130"]
        for dataSet in dataSetsToTest {
            theDataCentre.activeDataSetName = dataSet
            let (optimalWaveletPeakWidthKb, _) = theDataCentre.optimalWaveletWidthForSinglePeak(chromosome: chromosomeName, startKb: startPointKb, endKb: endPointKb)
            let waveletWidthBins = theDataCentre.waveletWidthBinsFromKb(waveletPeakWidthKb: optimalWaveletPeakWidthKb)
            let (_, waveletResult) = theDataCentre.analyseRegionWithWavelet(chromosome: chromosomeName, startKb: startPointKb, endKb: endPointKb, waveletPeakWidthKb: optimalWaveletPeakWidthKb)!
            let scaleFactorForDisplay = Double(waveletWidthBins) * 0.2                // scale the wavelet result so it is a bit below the replication signal
            let waveletResultScaled = waveletResult.map { $0 * scaleFactorForDisplay }
            let replicationPortion = Array(theDataCentre.activeDataSet.replicationSignalsDict[chromosomeName]!.repSignal[startBin...endBin])
            let waveletPeakHeightCutoff = minimumPeakHeightCutoffFromLateData(waveletPeakWidthKb: optimalWaveletPeakWidthKb, percentileCutoff: percentRemovalTextField.doubleValue) * scaleFactorForDisplay
            let cutOffLine: [(Double, Double)] = [(theDataCentre.binMidpointKbp(binNumber: startBin)/1000, waveletPeakHeightCutoff), (theDataCentre.binMidpointKbp(binNumber: endBin)/1000, waveletPeakHeightCutoff)]
            
            let resultGraph = JBLineGraph()
            resultGraph.createEquallySpacedYDataSetNamed("replication", withYData: replicationPortion, firstX: startPointKb/1000, XDelta: theDataCentre.dataBinSizeKb/1000)
            resultGraph.setPlotType(.barGraph, dataSetName: "replication")
            resultGraph.setPlotFillColor(NSColor(red: 0.4, green: 0.8, blue: 1, alpha: 1), dataSetName: "replication")
            resultGraph.setLineWidth(0, dataSetName: "replication")
            resultGraph.setBarLineWidth(0, dataSetName: "replication")
            resultGraph.createEquallySpacedYDataSetNamed("wavelet", withYData: waveletResultScaled, firstX: startPointKb/1000, XDelta: theDataCentre.dataBinSizeKb/1000)
            resultGraph.setPlotSymbolType(JBGraphSymbol.none, dataSetName: "wavelet")
            resultGraph.setLineColor(NSColor(red: 1, green: 0, blue: 0, alpha: 1), dataSetName: "wavelet")
            resultGraph.createXYDataSetNamed("cutoff", fromXYTupleArray: cutOffLine)
            resultGraph.setPlotSymbolType(JBGraphSymbol.none, dataSetName: "cutoff")
            resultGraph.setLineDashType(lineDashType.fineHalfDash, dataSetName: "cutoff")
            resultGraph.graphTitle = String(format: "optimal wavelet peak width %.0f kb cutoff %0.f%%\n%@", optimalWaveletPeakWidthKb, percentRemovalTextField.doubleValue, dataSetString)
            resultGraph.theXGraphAxis.userDefinedAxisVariables.axisMin = Decimal(startPointTextField.doubleValue)
            resultGraph.theXGraphAxis.userDefinedAxisVariables.axisMax = Decimal(endPointTextField.doubleValue)
            resultGraph.theYGraphAxis.userDefinedAxisVariables.axisMin = Decimal(0)
            resultGraph.setXAxisTitle("chromosome " + chromosomeName + " (Mbp)")
            resultGraph.setYAxisTitle("replication signal")
            resultGraph.createGraph(displayGraph: true)
            
            theDataCentre.activeDataSetName = dataSetPopupButton.selectedItem!.title
        }
    }
    
    @IBAction func singlePeakGaussFit(sender: NSButton) {
        if !checkChromosomeParameters() {
            return
        }
        let chromosomeName = chromosomeNameTextField.stringValue
        let startPointKb = startPointTextField.doubleValue * 1000
        let startBin = Int(startPointKb / theDataCentre.dataBinSizeKb)
        let endPointKb = endPointTextField.doubleValue * 1000
        let endBin = Int(endPointKb / theDataCentre.dataBinSizeKb)
        let dataSetsToTest = ["TM_10-40", "TM_40-70", "TM_70-100", "TM_100-130"]
        for dataSet in dataSetsToTest {
            theDataCentre.activeDataSetName = dataSet
            let replicationPortion = Array(theDataCentre.activeDataSet.replicationSignalsDict[chromosomeName]!.repSignal[startBin...endBin])
            let GaussParams = theDataCentre.optimalGaussianForSinglePeak(chromosome: chromosomeName, startKb: startPointKb, endKb: endPointKb)
            var GaussCurve = [Double]()
            let theGaussian = Gaussian(parameters: GaussParams)
            for index in 0..<replicationPortion.count {
                GaussCurve.append(theGaussian.valueFor(x: Double(index)))
            }
            let fullWidthAtHalfMax = theGaussian.fullWidthAtHalfMaximum() * theDataCentre.dataBinSizeKb
            let centrePointMb = (startPointKb + (theGaussian.centrePoint * theDataCentre.dataBinSizeKb)) / 1000
            
            let resultGraph = JBLineGraph()
            resultGraph.createEquallySpacedYDataSetNamed("replication", withYData: replicationPortion, firstX: startPointKb/1000, XDelta: theDataCentre.dataBinSizeKb/1000)
            resultGraph.setPlotType(.barGraph, dataSetName: "replication")
            resultGraph.setPlotFillColor(NSColor(red: 0.4, green: 0.8, blue: 1, alpha: 1), dataSetName: "replication")
            resultGraph.setLineWidth(0, dataSetName: "replication")
            resultGraph.setBarLineWidth(0, dataSetName: "replication")
            resultGraph.createEquallySpacedYDataSetNamed("fitted", withYData: GaussCurve, firstX: startPointKb/1000, XDelta: theDataCentre.dataBinSizeKb/1000)
            resultGraph.setPlotSymbolType(JBGraphSymbol.none, dataSetName: "fitted")
            resultGraph.theXGraphAxis.userDefinedAxisVariables.axisMin = Decimal(startPointTextField.doubleValue)
            resultGraph.theXGraphAxis.userDefinedAxisVariables.axisMax = Decimal(endPointTextField.doubleValue)
            resultGraph.theYGraphAxis.userDefinedAxisVariables.axisMin = Decimal(0)
            resultGraph.setXAxisTitle("chromosome " + chromosomeName + " (Mbp)")
            resultGraph.setYAxisTitle("replication signal")
            resultGraph.graphTitle = String(format: "Gaussian fit, %@\nheight %.0f, centre %.1f, stdev %.1f, FWHM %.0f", dataSetString, GaussParams[0], centrePointMb, GaussParams[2], fullWidthAtHalfMax)
            resultGraph.createGraph(displayGraph: true)
        }
        theDataCentre.activeDataSetName = dataSetPopupButton.selectedItem!.title
    }
    
    
   
    @IBAction func isolatedPeakSearch(sender: NSButton) {
        let waveletPeakWidthKb = waveletWidthTextField.doubleValue
        let (earlyTDs, _) = theDataCentre.U2OSTimings.allDomains()
        let waveletPeakHeightCutoff = minimumPeakHeightCutoffFromLateData(waveletPeakWidthKb: waveletPeakWidthKb, percentileCutoff: percentRemovalTextField.doubleValue)
        let isolatedPeaks = theDataCentre.isolatedPeaksInTimingDomains(timingDomainArray: earlyTDs, waveletPeakWidthKb: waveletPeakWidthKb, minimumHeight: waveletPeakHeightCutoff)
        var outputString = ""
        for aPeak in isolatedPeaks {
            outputString += String(format: "%@\t%.0f\n", aPeak.chromosome, aPeak.peakPositionKb)
        }
        
        let savePanel = NSSavePanel()
        savePanel.nameFieldStringValue = "isolated peaks list"
        savePanel.runModal()
        var panelResult = savePanel.url
        if panelResult == nil { return }
        panelResult = panelResult!.deletingPathExtension()
        panelResult = panelResult!.appendingPathExtension("txt")
        let dataPathString = panelResult!.path
        try? outputString.write(toFile: dataPathString, atomically: true, encoding: String.Encoding.utf8)
    }
    
    @IBAction func earlyVLateRepSignals(sender: NSButton) {
        let (earlyTDs, lateTDs) = theDataCentre.U2OSTimings.allDomains()
        let earlyRepSignals = theDataCentre.totalRepSignalPerTimingDomain(timingDomains: earlyTDs)
        let lateRepSignals = theDataCentre.totalRepSignalPerTimingDomain(timingDomains: lateTDs)
        let totalEarlyRep = earlyRepSignals.reduce(0, {$0 + $1}) / 1000000
        let totalLateRep = lateRepSignals.reduce(0, {$0 + $1}) / 1000000
        let earlyLateRepSignalsGraph = JBLineGraph()
        earlyLateRepSignalsGraph.createSequentialYDataSetNamed(nil, withYData: [totalEarlyRep, totalLateRep])
        earlyLateRepSignalsGraph.setPlotType(.barGraph, dataSetName: nil)
        earlyLateRepSignalsGraph.setXAxisTitle("early - late")
        earlyLateRepSignalsGraph.setYAxisTitle("total replication signal x 10^6")
        earlyLateRepSignalsGraph.graphTitle = String(format: "early signal %.1f x 10^6, late signal %.1f x 10^6\n%@",totalEarlyRep, totalLateRep, dataSetString)
        earlyLateRepSignalsGraph.createGraph(displayGraph: true)
    }
    
    @IBAction func peakWaveletWidthSweep(sender: NSButton) {
        let peakWidthsToTest: [Double] = [100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 1000, 1500, 3000]
        let (earlyTDs, _) = theDataCentre.U2OSTimings.allDomains()
        var medianPeakHeights = [Double]()
        var meanPeakHeights = [Double]()
        var peakHeightSums = [Double]()
        var peakNumber = [Double]()
        for waveletPeakWidthKb in peakWidthsToTest {
            let waveletPeakHeightCutoff = minimumPeakHeightCutoffFromLateData(waveletPeakWidthKb: waveletPeakWidthKb, percentileCutoff: percentRemovalTextField.doubleValue)
            var earlyPeakHeights = theDataCentre.waveletPeakHeightsInTimingDomains(timingDomainArray: earlyTDs, waveletPeakWidthKb: waveletPeakWidthKb)
            earlyPeakHeights = earlyPeakHeights.filter { $0 >= waveletPeakHeightCutoff }
            earlyPeakHeights.sort { $0 < $1 }
            let earlyPeakHeightMedian = earlyPeakHeights[earlyPeakHeights.count / 2]
            medianPeakHeights.append(earlyPeakHeightMedian)
            let earlyPeakHeightSum = earlyPeakHeights.reduce(0, { $0 + $1 })
            peakHeightSums.append(earlyPeakHeightSum)
            let earlyPeakHeightMean = earlyPeakHeightSum / Double(earlyPeakHeights.count)
            meanPeakHeights.append(earlyPeakHeightMean)
            peakNumber.append(Double(earlyPeakHeights.count))
        }
    
        let meanHeightsGraph = JBLineGraph()
        meanHeightsGraph.createXYDataSetNamed(nil, Xdata: peakWidthsToTest, YData: meanPeakHeights)
        meanHeightsGraph.setXAxisTitle("wavelet peak width")
        meanHeightsGraph.setYAxisTitle("mean peak height")
        meanHeightsGraph.graphTitle = String(format: "mean peak heights, cutoff %.0f%%\n%@", percentRemovalTextField.doubleValue, dataSetString)
        meanHeightsGraph.createGraph(displayGraph: true)
        
        let medianHeightsGraph = JBLineGraph()
        medianHeightsGraph.createXYDataSetNamed(nil, Xdata: peakWidthsToTest, YData: medianPeakHeights)
        medianHeightsGraph.setXAxisTitle("wavelet peak width")
        medianHeightsGraph.setYAxisTitle("median peak height")
        medianHeightsGraph.graphTitle = String(format: "median peak heights, cutoff %.0f%%\n%@", percentRemovalTextField.doubleValue, dataSetString)
        medianHeightsGraph.createGraph(displayGraph: true)
        
        let heightSumsGraph = JBLineGraph()
        heightSumsGraph.createXYDataSetNamed(nil, Xdata: peakWidthsToTest, YData: peakHeightSums)
        heightSumsGraph.setXAxisTitle("wavelet peak width")
        heightSumsGraph.setYAxisTitle("peak sum")
        heightSumsGraph.graphTitle = String(format: "peak height sums, cutoff %.0f%%\n%@", percentRemovalTextField.doubleValue, dataSetString)
        heightSumsGraph.createGraph(displayGraph: true)
        
        let peakNumberGraph = JBLineGraph()
        peakNumberGraph.createXYDataSetNamed(nil, Xdata: peakWidthsToTest, YData: peakNumber)
        peakNumberGraph.setXAxisTitle("wavelet peak width")
        peakNumberGraph.setYAxisTitle("number of peaks")
        peakNumberGraph.graphTitle = String(format: "peak number, cutoff %.0f%%\n%@", percentRemovalTextField.doubleValue, dataSetString)
        peakNumberGraph.createGraph(displayGraph: true)
    }
    
    @IBAction func peakPercentileHeightSweep(sender: NSButton) {
        let waveletPeakWidthKb = waveletWidthTextField.doubleValue
        let (earlyTDs, _) = theDataCentre.U2OSTimings.allDomains()
        let waveletPeaksInEarlyDNA = theDataCentre.waveletPeaksInTimingDomains(timingDomainArray: earlyTDs, waveletPeakWidthKb: waveletPeakWidthKb)
        var earlyPeakHeights = waveletPeaksInEarlyDNA.map {$0.peakHeight}
        earlyPeakHeights.sort { $0 < $1}
        var medianHeights = [Double]()
        var meanHeights = [Double]()
        var peaks = [Double]()
        
        for percentile in 0...100 {
            let waveletPeakHeightCutoff = minimumPeakHeightCutoffFromLateData(waveletPeakWidthKb: waveletPeakWidthKb, percentileCutoff: Double(percentile))
            let thePeakHeights = earlyPeakHeights.filter { $0 > waveletPeakHeightCutoff }
            medianHeights.append(thePeakHeights[thePeakHeights.count/2])
            let peakHeightSum = thePeakHeights.reduce(0, { $0 + $1 })
            meanHeights.append(peakHeightSum / Double(thePeakHeights.count))
            peaks.append(Double(thePeakHeights.count))
        }
        
        let medianHeightsGraph = JBLineGraph()
        medianHeightsGraph.createEquallySpacedYDataSetNamed(nil, withYData: medianHeights, firstX: 0, XDelta: 1)
        medianHeightsGraph.setPlotSymbolType(.none, dataSetName: nil)
        medianHeightsGraph.setXAxisTitle("percentile cutoff")
        medianHeightsGraph.setYAxisTitle("median peak height")
        medianHeightsGraph.graphTitle = String(format: "median peak heights, wavelet width %.0f kb\n%@", waveletPeakWidthKb, dataSetString)
        medianHeightsGraph.createGraph(displayGraph: true)
        
        let meanHeightsGraph = JBLineGraph()
        meanHeightsGraph.createEquallySpacedYDataSetNamed(nil, withYData: meanHeights, firstX: 0, XDelta: 1)
        meanHeightsGraph.setPlotSymbolType(.none, dataSetName: nil)
        meanHeightsGraph.setXAxisTitle("percentile cutoff")
        meanHeightsGraph.setYAxisTitle("mean peak height")
        meanHeightsGraph.graphTitle = String(format: "mean peak heights, wavelet width %.0f kb\n%@", waveletPeakWidthKb, dataSetString)
        meanHeightsGraph.createGraph(displayGraph: true)
        
        let peakNumberGraph = JBLineGraph()
        peakNumberGraph.createEquallySpacedYDataSetNamed(nil, withYData: peaks, firstX: 0, XDelta: 1)
        peakNumberGraph.setPlotSymbolType(.none, dataSetName: nil)
        peakNumberGraph.setXAxisTitle("percentile cutoff")
        peakNumberGraph.setYAxisTitle("number of peaks")
        peakNumberGraph.graphTitle = String(format: "peak number, wavelet width %.0f kb\n%@", waveletPeakWidthKb, dataSetString)
        peakNumberGraph.createGraph(displayGraph: true)
    }

    @IBAction func earlyDomainPeakHeights(sender: NSButton) {
        guard theDataCentre.activeDataSet.replicationSignalsDict.count > 0 else {
            return
        }
        let (earlyTDs, _) = theDataCentre.U2OSTimings.allDomains()
        let waveletPeakWidthKb = waveletWidthTextField.doubleValue
        var earlyPeakHeights = theDataCentre.waveletPeakHeightsInTimingDomains(timingDomainArray: earlyTDs, waveletPeakWidthKb: waveletPeakWidthKb)
        let waveletPeakHeightCutoff = minimumPeakHeightCutoffFromLateData(waveletPeakWidthKb: waveletPeakWidthKb, percentileCutoff: percentRemovalTextField.doubleValue)
        earlyPeakHeights = earlyPeakHeights.filter { $0 >= waveletPeakHeightCutoff }
        earlyPeakHeights.sort { $0 < $1 }
        let earlyPeakHeightMedian = earlyPeakHeights[earlyPeakHeights.count / 2]
        let earlyPeakHeightMean = earlyPeakHeights.reduce(0, { $0 + $1 }) / Double(earlyPeakHeights.count)
        let earlyHeightsGraph = JBLineGraph()
        earlyHeightsGraph.createFrequencyDataSetNamed(nil, Data: earlyPeakHeights, binSize: 100)
        earlyHeightsGraph.setXAxisTitle("peak height")
        earlyHeightsGraph.setYAxisTitle("frequency")
        earlyHeightsGraph.graphTitle = String(format: "%d early peaks, height mean %.0f median %.0f \nwavelet peak width %.0f, cutoff %.0f%% at %.0f\n%@", earlyPeakHeights.count, earlyPeakHeightMean, earlyPeakHeightMedian, waveletPeakWidthKb, percentRemovalTextField.doubleValue, waveletPeakHeightCutoff, dataSetString)
        earlyHeightsGraph.createGraph(displayGraph: true)
    }
        
    @IBAction func waveletPeakSeparations(sender: NSButton) {
        let waveletPeakWidthKb = waveletWidthTextField.doubleValue
        let (earlyTDs, _) = theDataCentre.U2OSTimings.allDomains()
        let waveletPeakHeightCutoff = minimumPeakHeightCutoffFromLateData(waveletPeakWidthKb: waveletPeakWidthKb, percentileCutoff: percentRemovalTextField.doubleValue)
        var earlyPeakSeparations = theDataCentre.peakSeparations(waveletPeakWidthKb: waveletWidthTextField.doubleValue, minimumPeakHeight: waveletPeakHeightCutoff, inTimingZones: earlyTDs)
        earlyPeakSeparations.sort {$0 < $1}
        let earlySeparationsMedian = earlyPeakSeparations[earlyPeakSeparations.count/2]
        let earlySeparationsMean = earlyPeakSeparations.reduce(0, { $0 + $1 }) / Double(earlyPeakSeparations.count)
        let earlySeparationsVariance = earlyPeakSeparations.reduce(0, { $0 + ($1-earlySeparationsMean)*($1-earlySeparationsMean) }) / Double(earlyPeakSeparations.count)
        let earlySeparationsStdv = sqrt(earlySeparationsVariance)
        let earlySeparationsGraph = JBLineGraph()
        earlySeparationsGraph.createFrequencyDataSetNamed(nil, Data: earlyPeakSeparations, binSize: 100)
        earlySeparationsGraph.setXAxisTitle("peak separation (kb)")
        earlySeparationsGraph.setYAxisTitle("frequency")
        earlySeparationsGraph.graphTitle = String(format: "%d peak seps, mean %.0f, median %.0f, stdev %.1f\nwavelet width %.0f kb, cutoff %.0f%% \n%@ ()", earlyPeakSeparations.count, earlySeparationsMean, earlySeparationsMedian, earlySeparationsStdv, waveletPeakWidthKb, percentRemovalTextField.doubleValue, dataSetString)
        earlySeparationsGraph.createGraph(displayGraph: true)
    }

    @IBAction func peaksPerEarlyTimingDomain(sender: NSButton) {
        // returns the number of peaks found in each timing domain in a tuble with the size of the domain
        let waveletPeakWidthKb = waveletWidthTextField.doubleValue
        let (earlyTDs, _) = theDataCentre.U2OSTimings.allDomains()
        let waveletPeakHeightCutoff = minimumPeakHeightCutoffFromLateData(waveletPeakWidthKb: waveletPeakWidthKb, percentileCutoff: percentRemovalTextField.doubleValue)
        let peaksPerTD = theDataCentre.peaksPerTimingDomain(waveletPeakWidthKb: waveletPeakWidthKb, minimumPeakHeight: waveletPeakHeightCutoff, inTimingDomains: earlyTDs)
        let peaksPerTDAsDoubles = peaksPerTD.map { ($0.sizeOfTimingDomainKb/1000, Double($0.numberOfPeaks)) }
        let peaksPerTDGraph = JBLineGraph()
        peaksPerTDGraph.createXYDataSetNamed(nil, fromXYTupleArray: peaksPerTDAsDoubles)
        peaksPerTDGraph.setLineWidth(0, dataSetName: nil)
        peaksPerTDGraph.setXAxisTitle("size of early timing domain (Mbp)")
        peaksPerTDGraph.setYAxisTitle("number of peaks")
        peaksPerTDGraph.graphTitle = String(format: "peaks per timing domain \nwavelet width %.0f kb, cutoff %.0f%% \n%@ ()", waveletPeakWidthKb, percentRemovalTextField.doubleValue, dataSetString)
        peaksPerTDGraph.createGraph(displayGraph: true)
    }
    
    @IBAction func valleyFilling(sender: NSButton) {
        // identifies 'valleys' in the TM_10-40 data set, then iterates thoguh all the valleys recording their mean replication signal intensity at all four timepoints
        let edgeOffsetDueToForkElongationKb = 120 * 1.5        // the amount the edge of each peak is extended by fork movement over 2 hours
        let peakPostionError = 100.0                              // peak position can shift ± 100 kb and still be the 'same' peak
        let waveletPeakWidthKb = waveletWidthTextField.doubleValue
        theDataCentre.activeDataSetName = "TM_10-40"
        let firstTimepointValleyList = theDataCentre.allValleys(waveletPeakWidthKb: waveletPeakWidthKb, edgeOffsetDueToForkElongationKb: edgeOffsetDueToForkElongationKb)
        let dataSetsToTest = ["TM_10-40", "TM_40-70", "TM_70-100", "TM_100-130"]
        var consistentValleyPercentageDensities = [String : [Double]]()   // valleys that exist in all 4 timepoints; the densities are the means densities of the valleys as a percent of flanking peaks
        var consistentValleyMinimaAsPercent = [String : [Double]]()   // valleys that exist in all 4 timepoints; this stores the replication signal minima as a percent of flanking peaks
        for aDataSet in dataSetsToTest {
            consistentValleyPercentageDensities[aDataSet] = [Double]()
            consistentValleyMinimaAsPercent[aDataSet] = [Double]()
        }
        for aValley in firstTimepointValleyList {             // remove valleys where the flanking peaks do not persist across all timepoints, then get their fractional intensities
            let theChromosome = aValley.leftPeak.chromosome
            let leftPeakArray = theDataCentre.allPeaksAtPosition(chromosome: theChromosome, peakPositionKb: aValley.leftPeak.peakPositionKb, positionErrorKb: peakPostionError, inDataSets: dataSetsToTest, waveletPeakWidthKb: waveletPeakWidthKb)
            if leftPeakArray.reduce(0, { $1 != nil ? $0 + 1 : $0 }) < dataSetsToTest.count { continue }    // lose it if the left peak isn't always there
            let rightPeakArray = theDataCentre.allPeaksAtPosition(chromosome: theChromosome, peakPositionKb: aValley.rightPeak.peakPositionKb, positionErrorKb: peakPostionError, inDataSets: dataSetsToTest, waveletPeakWidthKb: waveletPeakWidthKb)
            if rightPeakArray.reduce(0, { $1 != nil ? $0 + 1 : $0 }) < dataSetsToTest.count { continue }        // lose it if the right peak isn't always there
            let valleyLeftEdgeKb = aValley.leftPeak.endKb + edgeOffsetDueToForkElongationKb                     // valley edges defined by the first time point
            let valleyRightEdgeKb = aValley.rightPeak.startKb - edgeOffsetDueToForkElongationKb
            for timepointIndex in 0..<dataSetsToTest.count {
                let aDataSet = dataSetsToTest[timepointIndex]
                theDataCentre.activeDataSetName = aDataSet
                let (_, _, minimumValleyIntensity, meanValleyIntensity) = theDataCentre.repSignalMetricsInRegion(chromosome: theChromosome, startKb: valleyLeftEdgeKb, endKb: valleyRightEdgeKb)
                let meanPeakIntensity = (leftPeakArray[timepointIndex]!.maxReplicationSignal + rightPeakArray[timepointIndex]!.maxReplicationSignal) / 2
                consistentValleyPercentageDensities[aDataSet]!.append(meanValleyIntensity * 100 / meanPeakIntensity)        // 100 because it's percent
                consistentValleyMinimaAsPercent[aDataSet]!.append(minimumValleyIntensity * 100 / meanPeakIntensity)        // 100 because it's percent)
            }
        }
        
        for aDataSet in dataSetsToTest {                // plot graphs for mean and minimum valley percentages
            let meanRepSignalMean = consistentValleyPercentageDensities[aDataSet]!.reduce(0, {$0 + $1}) / Double(consistentValleyPercentageDensities[aDataSet]!.count)
            let meanRepSignalVariance = consistentValleyPercentageDensities[aDataSet]!.reduce(0, { $0 + ($1-meanRepSignalMean)*($1-meanRepSignalMean) }) / Double(consistentValleyPercentageDensities[aDataSet]!.count)
            let meanRepSignalStdev = sqrt(meanRepSignalVariance)
            let valleyMeansGraph = JBLineGraph()
            valleyMeansGraph.createFrequencyDataSetNamed(nil, Data: consistentValleyPercentageDensities[aDataSet]!, binSize: 10)
            valleyMeansGraph.setXAxisTitle("valley mean signal (% flank)")
            valleyMeansGraph.setYAxisTitle("frequency")
            valleyMeansGraph.setPlotType(.lineGraph, dataSetName: nil)
            valleyMeansGraph.graphTitle = String(format: "%d valley means (%d at 1st time), mean %.4f\n stdev %.4f, wavelet width %.0f kb, \n %@ %@ ()", consistentValleyPercentageDensities[aDataSet]!.count, firstTimepointValleyList.count, meanRepSignalMean, meanRepSignalStdev, waveletPeakWidthKb, dataFileName!, aDataSet)
            valleyMeansGraph.createGraph(displayGraph: true)
            
            let minRepSignalMean = consistentValleyMinimaAsPercent[aDataSet]!.reduce(0, {$0 + $1}) / Double(consistentValleyMinimaAsPercent[aDataSet]!.count)
            let minRepSignalVariance = consistentValleyMinimaAsPercent[aDataSet]!.reduce(0, { $0 + ($1-minRepSignalMean)*($1-minRepSignalMean) }) / Double(consistentValleyMinimaAsPercent[aDataSet]!.count)
            let minRepSignalStdev = sqrt(minRepSignalVariance)
            let valleyMinimaGraph = JBLineGraph()
            valleyMinimaGraph.createFrequencyDataSetNamed(nil, Data: consistentValleyMinimaAsPercent[aDataSet]!, binSize: 10)
            valleyMinimaGraph.setXAxisTitle("valley minimum signal (% flank)")
            valleyMinimaGraph.setYAxisTitle("frequency")
            valleyMinimaGraph.setPlotType(.lineGraph, dataSetName: nil)
            valleyMinimaGraph.graphTitle = String(format: "%d valley minima (%d at 1st time), mean %.4f\n stdev %.4f, wavelet width %.0f kb, \n %@ %@ ()", consistentValleyMinimaAsPercent[aDataSet]!.count, firstTimepointValleyList.count, minRepSignalMean, minRepSignalStdev, waveletPeakWidthKb, dataFileName!, aDataSet)
            valleyMinimaGraph.createGraph(displayGraph: true)
        }
        theDataCentre.activeDataSetName = dataSetPopupButton.selectedItem!.title
    }
    
    @IBAction func peakActivationSequence(sender: NSButton) {
        let maximumDistanceBetweenPeaksKb = 1600.0
        var sizeRankingsOfAdjacentPeaks = theDataCentre.sizeSequenceOfAdjacentPeaks(waveletPeakWidthKb: waveletWidthTextField.doubleValue, maximumDistanceBetweenPeaksKb: maximumDistanceBetweenPeaksKb)
        sizeRankingsOfAdjacentPeaks.sort(by: { $0.count < $1.count })
        
        let groupsOfThree = sizeRankingsOfAdjacentPeaks.filter({ $0.count == 3 })
        let tripPossibilityMatrix = [[1,2,3],[1,3,2],[2,1,3],[3,1,2]]
        var tripResultArray = Array(repeating: 0, count: tripPossibilityMatrix.count)
        for anArray in groupsOfThree {
            for index in 0..<tripPossibilityMatrix.count {
                if anArray == tripPossibilityMatrix[index] {
                    tripResultArray[index] += 1
                    break
                }
            }
        }
        
        var outputString = String(format: "peak activation order for wavelet width %.0f kb, cutoff %.0f%% \n%@\nmaximum separation between peaks %.0f kb\n\n", waveletWidthTextField.doubleValue, percentRemovalTextField.doubleValue, dataSetString, maximumDistanceBetweenPeaksKb)
        outputString += String(format: "%d triplets\n", groupsOfThree.count)
        for index in 0..<tripResultArray.count {
            outputString += String(format: "%d * %d-%d-%d\n", tripResultArray[index], tripPossibilityMatrix[index][0],tripPossibilityMatrix[index][1],tripPossibilityMatrix[index][2])
        }
        
        let groupsOfFour = sizeRankingsOfAdjacentPeaks.filter({ $0.count == 4 })
        let quadPossibilityMatrix = [[1,2,3,4],[1,2,4,3],[1,3,2,4],[1,3,4,2],[1,4,2,3],[1,4,3,2],[2,1,3,4],[2,1,4,3],[3,1,2,4],[3,1,4,2],[4,1,2,3],[4,1,3,2]]
        var quadResultArray = Array(repeating: 0, count: quadPossibilityMatrix.count)
        for anArray in groupsOfFour {
            for index in 0..<quadPossibilityMatrix.count {
                if anArray == quadPossibilityMatrix[index] {
                    quadResultArray[index] += 1
                    break
                }
            }
        }
        
        outputString += String(format: "\n%d quads\n", groupsOfFour.count)
        for index in 0..<quadResultArray.count {
            outputString += String(format: "%d * %d-%d-%d-%d\n", quadResultArray[index], quadPossibilityMatrix[index][0],quadPossibilityMatrix[index][1],quadPossibilityMatrix[index][2],quadPossibilityMatrix[index][3])
        }
        
        let maximumGroupSize = sizeRankingsOfAdjacentPeaks.reduce(0, { $1.count > $0 ? $1.count : $0 } )
        outputString += "\ngroup size frequncy\n"
        for index in 1...maximumGroupSize {
            let frequency = sizeRankingsOfAdjacentPeaks.reduce(0, { $1.count == index ? $0 + 1 : $0 })
            outputString += String(format: "%d: %d\n", index, frequency)
        }
        let savePanel = NSSavePanel()
        savePanel.nameFieldStringValue = "peak order"
        savePanel.runModal()
        var panelResult = savePanel.url
        if panelResult == nil { return }
        panelResult = panelResult!.deletingPathExtension()
        panelResult = panelResult!.appendingPathExtension("txt")
        let dataPathString = panelResult!.path
        try? outputString.write(toFile: dataPathString, atomically: true, encoding: String.Encoding.utf8)
    }
    
    @IBAction func peakNeighbourSimilarity(sender: NSButton) {
        let maximumDistanceBetweenPeaksKb = 1600.0
        let adjacentWaveletPeaks = theDataCentre.adjacentWaveletPeaks(waveletPeakWidthKb: waveletWidthTextField.doubleValue, maximumDistanceBetweenPeaksKb: maximumDistanceBetweenPeaksKb)
        let adjacentWaveletPeaksFiltered = adjacentWaveletPeaks.filter( { $0.count > 3 })       // at 3 or below, no difference between random and anti-ordered
        
        var degreeOfOrdering = theDataCentre.similarityInRepSignalsBetweenAdjacentPeaks(groupOfAdjacentWavelets: adjacentWaveletPeaksFiltered)
        degreeOfOrdering = degreeOfOrdering.map( {$0 >= 1 ? 0.99999 : $0 })          // stops silly rounding effects
        degreeOfOrdering = degreeOfOrdering.map( {$0 <= -1 ? -0.99999 : $0 })          // stops silly rounding effects
        let meanDegreeOfOrdering = degreeOfOrdering.reduce(0, {$0 + $1}) / Double(degreeOfOrdering.count)
        let degreeOfOrderingSorted = degreeOfOrdering.sorted(by: {$0 < $1})
        let medianDegreeOfOrdering = degreeOfOrderingSorted[degreeOfOrderingSorted.count/2]
        
        let theDeviationGraph = JBLineGraph()
        theDeviationGraph.createFrequencyDataSetNamed(nil, Data: degreeOfOrdering, binSize: 0.2)
        theDeviationGraph.setXAxisTitle("degree of peak ordering")
        theDeviationGraph.setYAxisTitle("frequency")
        theDeviationGraph.graphTitle = String(format: "%d peak separations, mean %.3f, median %.3f\nwavelet width %.0f kb, cutoff %.0f%% \n%@ ()", degreeOfOrdering.count, meanDegreeOfOrdering, medianDegreeOfOrdering, waveletWidthTextField.doubleValue, percentRemovalTextField.doubleValue, dataSetString)
        theDeviationGraph.createGraph(displayGraph: true)

    }
    
}
