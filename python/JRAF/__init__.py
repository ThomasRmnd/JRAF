import Sniper as sn
import types
sn.loadDll("libJRAF.so")
sn.loadDll("libCLHEPDict.so")

def useEventBuilder(self, name):
    eventbuilder = self.createTool(name)
    self.property("EventBuilder").set(name)
    self.eventbuilder = eventbuilder

def useMuonTagger(self, name):
    muontagger = self.createTool(name)
    self.eventbuilder.property("MuonTagger").set(name)
    self.muontagger = muontagger

def useLoader(self, name):
    loader = self.createTool(name)
    self.eventbuilder.property("Loader").set(name)
    self.loader = loader

def useCdFiller(self, name):
    cdfiller = self.createTool(name)
    self.eventbuilder.property("CdFiller").set(name)
    self.cdfiller = cdfiller

def useWpFiller(self, name):
    wpfiller = self.createTool(name)
    self.eventbuilder.property("WpFiller").set(name)
    self.wpfiller = wpfiller

def useTtFiller(self, name):
    ttfiller = self.createTool(name)
    self.eventbuilder.property("TtFiller").set(name)
    self.ttfiller = ttfiller

def useRecTool(self, name):
    rectool = self.createTool(name)
    self.eventbuilder.property("RecTool").set(name)
    self.rectool = rectool

def createAlg(task, name="JRAF"):
    alg = task.createAlg(name)
    alg.useEventBuilder = types.MethodType(useEventBuilder, alg)
    alg.useMuonTagger = types.MethodType(useMuonTagger, alg)
    alg.useLoader = types.MethodType(useLoader, alg)
    alg.useCdFiller = types.MethodType(useCdFiller, alg)
    alg.useWpFiller = types.MethodType(useWpFiller, alg)
    alg.useTtFiller = types.MethodType(useTtFiller, alg)
    alg.useRecTool = types.MethodType(useRecTool, alg)
    return alg