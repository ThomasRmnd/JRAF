import Sniper as sn
import types
sn.loadDll("libAnalysisGroupC.so")
sn.loadDll("libCLHEPDict.so")

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
    alg.useRecTool = types.MethodType(useRecTool, alg)
    alg.useClassifyTool = types.MethodType(useClassifyTool, alg)
    return alg