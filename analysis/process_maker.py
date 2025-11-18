import csv


input_txt = "run3datasetALL.txt"
output_py = "process.py"
xsec_csv  = "run3_datasets_XS.csv"

processes = {}
########################
## Read run3 xs file
########################
xsec_dict = {}
with open(xsec_csv, 'r') as f:
    reader = csv.reader(f)
    next(reader)  # header skip
    for row in reader:
        if len(row) < 2:
            continue
        name, xsec = row[0].strip(), row[1].strip()
        try:
            xsec = float(xsec)
        except:
            xsec = "XXX"
        xsec_dict[name] = xsec
        print('csv name', name)


########################
## Read run3 all dataset file
########################

with open(input_txt, 'r') as f:
    for line in f:
         line = line.strip()
         process_name = line # It will dataset name.
         base_name = process_name.split('_Tune')[0]
         short_tag = process_name.split('_')[0]
         print('process_name', process_name)

         processes[process_name] = {
             'base': base_name,
             'tag': short_tag,
             'xs' : xsec_dict[process_name]
         }

with open(output_py, 'w') as f:
    f.write("#!/usr/bin/env python\n\n")
    f.write("processes = {\n")
    f.write("\n    ## Data\n")
    f.write("    'JetMET':('JetMET','Data',1),\n")
    f.write("    'EGamma':('EGamma','Data',1),\n")
    current_tag = None
    for key, info in processes.items():
        if info['tag'] != current_tag:
            f.write(f"\n    ## {info['tag']}\n")
            current_tag = info['tag']
        f.write(f"    '{key}': ('{info['base']}', 'MC', {float(info['xs'])}),\n")
    f.write("}\n")

print(f"✅ Saved {len(processes)} processes to {output_py}")

savelines = '''
with open(input_txt, 'r') as f:
     for line in f:
         line = line.strip()
         if not line or not line.startswith('/'):
             continue
         parts = line.split('/')
         process_name = parts[1] # It will dataset name.
         base_name = process_name.split('_Tune')[0]
         short_tag = process_name.split('_')[0]

         processes[process_name] = {
             'base': base_name,
             'tag': short_tag
         }
'''
