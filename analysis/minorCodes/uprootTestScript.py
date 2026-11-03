import uproot

f = uproot.open("/data/mc/privatemc/2023BPix/sig_Mphi-1000_Mchi-1000/NanoAODv12_private/260211_060055/0000/nano_1.root")
print("keys =>", f.keys())
print("f['Events'].num_entries => ", f["Events"].num_entries)

events = f["Events"]
print("events.keys() => ", events.keys())

MET_etr = events.arrays(["MET_pt"], entry_stop=10)
print("MET_pt 10 ent => ",MET_etr)
