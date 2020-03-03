import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        ['/store/mc/RunIIFall17MiniAODv2/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/100000/86DF496F-7AC2-E811-B3C5-0025905C2CD2.root','/store/mc/RunIIFall17MiniAODv2/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/100000/E25EEFCD-B0C0-E811-9072-AC1F6B0DE228.root','/store/mc/RunIIFall17MiniAODv2/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/100000/E8782068-86C3-E811-A178-001E67E6F855.root','/store/mc/RunIIFall17MiniAODv2/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/110000/6EE94B14-84D4-E811-BA56-0CC47A6C138A.root','/store/mc/RunIIFall17MiniAODv2/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/100000/067E786D-ABC3-E811-8455-1866DA85DC8B.root','/store/mc/RunIIFall17MiniAODv2/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/100000/26ED3F0D-D3C3-E811-BA27-0CC47AFC3D34.root','/store/mc/RunIIFall17MiniAODv2/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/100000/2E9F5B07-1CC4-E811-B508-0CC47A1DF806.root','/store/mc/RunIIFall17MiniAODv2/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/100000/320FEF3C-E3C4-E811-8ED4-B499BAABD28C.root']
    )
)

runOnMC  = True


process.demo = cms.EDAnalyzer('hlt',
                                    RunOnMC  = cms.bool(runOnMC),
                                    genSrc =  cms.InputTag("prunedGenParticles"),

)


process.p = cms.Path(process.demo)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("RStreeEDBR_pickup.root")
                                   )

