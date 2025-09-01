import Sniper as sn
import types
sn.loadDll("libAnalysisGroupC.so")
sn.loadDll("libCLHEPDict.so")

def useRecTool(self, name):
    rectool = self.createTool(name)
    self.property("RecTool").set(name)
    self.rectool = rectool

def createAlg(task, name="AnalysisGroupC"):
    alg = task.createAlg(name)
    alg.useRecTool = types.MethodType(useRecTool, alg)
    return alg