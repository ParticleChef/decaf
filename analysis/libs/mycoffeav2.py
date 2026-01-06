import awkward
from coffea.nanoevents.methods import base, vector, candidate, nanoaod
from coffea.nanoevents import NanoAODSchema, BaseSchema

#behavior = {}
#behavior.update(base.behavior)
#behavior.update(candidate.behavior)
#print(candidate.behavior)
#print(NanoAODSchema.behavior)
#my_behavior = dict(NanoAODSchema.behavior)
my_behavior = {}
my_behavior.update(nanoaod.behavior)

def _set_repr_name(classname):
    def namefcn(self):
        return classname

    # behavior[("__typestr__", classname)] = classname[0].lower() + classname[1:]
    my_behavior[classname].__repr__ = namefcn

@awkward.mixin_class(my_behavior)
class AK15SubJet(vector.PtEtaPhiMLorentzVector, base.NanoCollection, base.Systematic):
    """NanoAOD AK15 subjet object"""

    @property
    def matched_gen(self):
        return self._events().SubGenJetAK15._apply_global_index(self.genJetIdxG)

_set_repr_name("AK15SubJet")

@awkward.mixin_class(my_behavior)
class AK15Jet(vector.PtEtaPhiMLorentzVector, base.NanoCollection, base.Systematic):
    """NanoAOD large radius jet object"""

    LOOSE = 0
    "jetId bit position"
    TIGHT = 1
    "jetId bit position"
    TIGHTLEPVETO = 2
    "jetId bit position"

    @property
    def isLoose(self):
        """Returns a boolean array marking loose jets according to jetId index"""
        return (self.jetId & (1 << self.LOOSE)) != 0

    @property
    def isTight(self):
        """Returns a boolean array marking tight jets according to jetId index"""
        return (self.jetId & (1 << self.TIGHT)) != 0

    @property
    def isTightLeptonVeto(self):
        """Returns a boolean array marking tight jets with explicit lepton veto according to jetId index"""
        return (self.jetId & (1 << self.TIGHTLEPVETO)) != 0

    @property
    def subjets(self):
        return self._events().AK15PuppiSubJet._apply_global_index(self.subJetIdxG)

    @property
    def matched_gen(self):
        return self._events().GenJetAK15._apply_global_index(self.genJetIdxG)


_set_repr_name("AK15Jet")

class CustomNanoAODSchema(NanoAODSchema):
    mixins = {
        **NanoAODSchema.mixins,
        "AK15PuppiJet": "AK15Jet",
        "AK15PuppiSubJet": "AK15SubJet",
    }
    all_cross_references = {
        **NanoAODSchema.all_cross_references,
        "AK15PuppiJet_genJetIdx":  "GenJetAK15",
        "AK15PuppiJet_subJetIdx1": "AK15PuppiSubJet",  
        "AK15PuppiJet_subJetIdx2": "AK15PuppiSubJet",
        "AK15PuppiSubJet_genJetIdx": "SubGenJetAK15",
    }
    nested_items = {
        **NanoAODSchema.nested_items,
        "AK15PuppiJet_subJetIdxG": ["AK15PuppiJet_subJetIdx1G", "AK15PuppiJet_subJetIdx2G"]
    }
    def __init__(self, base_form):
        for key in base_form["contents"].copy():
            if 'AK15Puppi_' in key:
                popped = base_form["contents"].pop(key)
                base_form["contents"][key.replace('AK15Puppi_','AK15PuppiJet_')] = popped
            #if '_Jet' in key:
            #    popped = base_form["contents"].pop(key)
            #    base_form["contents"][key.replace('_Jet','Jet')] = popped
            #    #base_form["contents"][key.replace('_Jet','Jet')]['form_key'] = \
            #    #                    base_form["contents"][key.replace('_Jet','Jet')]['form_key'].replace('_Jet','Jet')
            #if '_Subjet' in key:
            #    popped = base_form["contents"].pop(key)
            #    base_form["contents"][key.replace('_Subjet','SubJet')] = popped
            #    #base_form["contents"][key.replace('_Subjet','SubJet')]['form_key'] = \
            #    #                    base_form["contents"][key.replace('_Subjet','SubJet')]['form_key'].replace('_Subjet','SubJet')
        super().__init__(base_form)

    @property
    def behavior(self):
        return my_behavior

