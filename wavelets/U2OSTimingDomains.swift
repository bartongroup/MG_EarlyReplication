//
//  U2OSTimingDomains.swift
//  ReplicationWavelet
//
//  Created by Julian Blow on 11/11/2021.
//

import Foundation
import Cocoa

enum TimingDomainType {
    case early, late, ambiguous
}

let maximumValidChromosomeName = 5

struct TimingDomainDescription {
    let chromosome: String
    let start: Int
    let end: Int
    let domainType: TimingDomainType
}

struct EarlyAndLateArrays {
    let earlyDomains : [TimingDomainDescription]
    let lateDomains : [TimingDomainDescription]
}

struct U2OSTimingDomains {
    let timingDomainsByChromosome : [String : EarlyAndLateArrays]
    let timingsByChromosome : [String : [Double]]
    let binSizeBp: Int
    
    init() {
        let openPanel = NSOpenPanel()
        openPanel.message = "U2OS timing domains file"
        openPanel.runModal()
        let panelResult = openPanel.url
        let dataPathString = panelResult!.path
        var (tempTimingDomains, tempTimings, binSize) = U2OSTimingDomains.readTimingDomainsFromFile(path: dataPathString)
        let chromNames = Array(tempTimingDomains.keys)
        for chromosomeName in chromNames {                        // remove the long names which are things we can't map to our simple genome
            if chromosomeName.count > maximumValidChromosomeName {
                tempTimingDomains[chromosomeName] = nil
                tempTimings[chromosomeName] = nil
            }
        }
        tempTimingDomains["chrY"] = nil                      // U2OS is female
        tempTimings["chrY"] = nil
        timingDomainsByChromosome = tempTimingDomains
        timingsByChromosome = tempTimings
        binSizeBp = binSize
    }
    
    func domainsMatching(domainType : TimingDomainType, minimumSizeKb: Int) -> [TimingDomainDescription] {
        // returns all the timing domains as specified by domainType and with a minimum size as specified
        let minimumSize = minimumSizeKb * 1000
        var timingDomainArray = [TimingDomainDescription]()
        for (_, timingArrays) in timingDomainsByChromosome {
            var targetArray = timingArrays.earlyDomains
            if domainType == .late {
                targetArray = timingArrays.lateDomains
            }
            for domain in targetArray {
                if (domain.end - domain.start) >= minimumSize {
                    timingDomainArray.append(domain)
                }
            }
        }
        return timingDomainArray
    }
    
    func allDomains() -> (earlyDomains: [TimingDomainDescription], lateDomains: [TimingDomainDescription]) {
        // returns arrays of early and late timing domains
        var earlyDomains = [TimingDomainDescription]()
        var lateDomains = [TimingDomainDescription]()
        for (_, arrays) in timingDomainsByChromosome {
            earlyDomains += arrays.earlyDomains
            lateDomains += arrays.lateDomains
        }
        return (earlyDomains, lateDomains)
    }
    
    func timingDomainLevelsInRegion(chromosome: String, startKb: Double, endKb: Double) -> [Double] {
        // returns an array of all timing domains in the specified region
        let startBin = Int(startKb * 1000) * binSizeBp
        let endBin = Int(endKb * 1000) / binSizeBp
        guard timingsByChromosome[chromosome] != nil else {
            return [Double]()
        }
        let timingsArray = timingsByChromosome[chromosome]!
        if endBin >= timingsArray.count {                           // likely have missing data at chromosome end
            let paddingArray = Array(repeating: 0.0, count: endBin + 1 - timingsArray.count)
            return Array(timingsArray[startBin...(timingsArray.count - 1)]) + paddingArray
        }
        return Array(timingsArray[startBin...endBin])
    }
    
    static func readTimingDomainsFromFile(path: String) -> ([String : EarlyAndLateArrays],  [String : [Double]], Int) {
        // reads timing domain data from a specified file. The format of the file is assumed to be correct and there is no format checking.
        let TDFile = try? String(contentsOfFile: path)
        var fileLines = TDFile!.components(separatedBy: "\n")
        fileLines = fileLines.filter { !$0.starts(with: "#") }
        fileLines = fileLines.filter { !$0.starts(with: " ") }
        fileLines = fileLines.filter { !$0.starts(with: "\r") }
        fileLines = fileLines.filter { !$0.isEmpty }
        var chromDomainsDict = [String : EarlyAndLateArrays]()
        var chromTimingsDict = [String : [Double]]()
        var currentChromosome = ""
        var earlyDomainArray = [TimingDomainDescription]()
        var lateDomainArray = [TimingDomainDescription]()
        var startOfCurrentDomain = 0
        var endOfCurrentDomain = 0
        var currentDomainType = TimingDomainType.ambiguous
        var newDomainType = TimingDomainType.ambiguous
        var timingsArray = [Double]()
        let firstLineComponents = fileLines[0].components(separatedBy: "\t")            // to extract bin size
        let binWidthBP = Int(firstLineComponents[2])! - Int(firstLineComponents[1])!
        var nextBinStart = 0
        var currentTimingValue = 0.0
        for line in fileLines {
            let lineComponents = line.components(separatedBy: "\t")
            if currentChromosome != lineComponents[0] {                 // found a new chromosome
                if currentDomainType == .early {                 // first close down the last domains of the chromosome - not elegant!!!
                    earlyDomainArray.append(TimingDomainDescription(chromosome: currentChromosome, start: startOfCurrentDomain, end: endOfCurrentDomain, domainType: .early))
                }
                else if currentDomainType == .late {
                    lateDomainArray.append(TimingDomainDescription(chromosome: currentChromosome, start: startOfCurrentDomain, end: endOfCurrentDomain, domainType: .late))
                }
                precondition(currentDomainType == .ambiguous || earlyDomainArray.count > 0 || lateDomainArray.count > 0, "finished chromosome with no domains")
                if earlyDomainArray.count > 0 || lateDomainArray.count > 0 {
                    precondition(chromDomainsDict[currentChromosome] == nil, "trying to add a chromosome that already exists")
                    chromDomainsDict[currentChromosome] = EarlyAndLateArrays(earlyDomains: earlyDomainArray, lateDomains: lateDomainArray)
                    chromTimingsDict[currentChromosome] = timingsArray
                }
                earlyDomainArray.removeAll()
                lateDomainArray.removeAll()
                timingsArray.removeAll()
                currentDomainType = TimingDomainType.ambiguous
                startOfCurrentDomain = Int(lineComponents[1])!
                currentChromosome = lineComponents[0]
                nextBinStart = 0
                currentTimingValue = 0.0
            }
            while nextBinStart < Int(lineComponents[1])! {
                timingsArray.append(currentTimingValue)                                  // fill in any gaps with current value
                nextBinStart += binWidthBP
            }
            nextBinStart += binWidthBP
            currentTimingValue = Double(lineComponents[3])!
            timingsArray.append(currentTimingValue)
            if currentTimingValue > 0 {
                newDomainType = .early
            } else {
                newDomainType = .late
            }
            if newDomainType != currentDomainType {
                if currentDomainType == .early {
                    precondition(endOfCurrentDomain > startOfCurrentDomain, "trying to create early domain with start > end")
                    earlyDomainArray.append(TimingDomainDescription(chromosome: currentChromosome, start: startOfCurrentDomain, end: endOfCurrentDomain, domainType: .early))
                }
                else if currentDomainType == .late {
                    precondition(endOfCurrentDomain > startOfCurrentDomain, "trying to create late domain with start > end")
                    lateDomainArray.append(TimingDomainDescription(chromosome: currentChromosome, start: startOfCurrentDomain, end: endOfCurrentDomain, domainType: .late))
                }
                startOfCurrentDomain = Int(lineComponents[1])!
                currentDomainType = newDomainType
            }
            endOfCurrentDomain = Int(lineComponents[2])!
            
        }
        if (earlyDomainArray.count > 0) || (lateDomainArray.count > 0) {                             // collect last domain (not elegant!!!!)
            precondition(chromDomainsDict[currentChromosome] == nil, "trying to add a chromosome that already exists")
            chromDomainsDict[currentChromosome] = EarlyAndLateArrays(earlyDomains: earlyDomainArray, lateDomains: lateDomainArray)
            chromTimingsDict[currentChromosome] = timingsArray
        }
        return (chromDomainsDict, chromTimingsDict, binWidthBP)
    }
    
}





