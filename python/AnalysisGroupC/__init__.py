import Sniper as sn
import types
sn.loadDll("libAnalysisGroupC.so")
sn.loadDll("libCLHEPDict.so")

def useLoader(self, name):
    loader = self.createTool(name)
    self.property("Loader").set(name)
    self.loader = loader

def useCdFiller(self, name):
    cdfiller = self.createTool(name)
    self.property("CdFiller").set(name)
    self.cdfiller = cdfiller

def useWpFiller(self, name):
    wpfiller = self.createTool(name)
    self.property("WpFiller").set(name)
    self.wpfiller = wpfiller

def useTtFiller(self, name):
    ttfiller = self.createTool(name)
    self.property("TtFiller").set(name)
    self.ttfiller = ttfiller

def useClassifyLoader(self):
    classify_loader = self.createTool("BasicLoader")
    classify_wpfiller = self.createTool("AkiraWpRangeFiller")
    self.classify_loader = classify_loader
    self.classify_wpfiller = classify_wpfiller

def useRecTool(self, name):
    rectool = self.createTool(name)
    self.property("RecTool").set(name)
    self.rectool = rectool

def useClassifyTool(self, name):
    classifytool = self.createTool(name)
    self.property("ClassifyTool").set(name)
    self.classifytool = classifytool

def createAlg(task, name="AnalysisGroupC"):
    alg = task.createAlg(name)
    alg.useLoader = types.MethodType(useLoader, alg)
    alg.useCdFiller = types.MethodType(useCdFiller, alg)
    alg.useWpFiller = types.MethodType(useWpFiller, alg)
    alg.useTtFiller = types.MethodType(useTtFiller, alg)
    alg.useClassifyLoader = types.MethodType(useClassifyLoader, alg)
    alg.useRecTool = types.MethodType(useRecTool, alg)
    alg.useClassifyTool = types.MethodType(useClassifyTool, alg)
    return alg