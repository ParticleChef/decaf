import pandas as pd
import glob

def compare_columns(file1):
    df1 = pd.read_csv(file1, delim_whitespace=True, names=["event number", "run number", "lumi number"], skiprows=1)
    df1.columns = [col.strip() for col in df1.columns]

    file_list = glob.glob("egmNumberlist/egm_gcr_EGamma*")
    print(file_list)
    all_data = []

    for files in file_list:
        split_df = pd.read_csv(files, delim_whitespace=True, names=["event number", "run number", "lumi number"], skiprows=1)
        split_df.columns = [col.strip() for col in split_df.columns]
        all_data.append(split_df)

    df2 = pd.concat(all_data, ignore_index=True).drop_duplicates()

    merged_df = pd.merge(df1, df2, how="outer", indicator=True)

    diff = merged_df[merged_df["_merge"] != "both"]

    with open("egm_compare.txt", "w") as out:
        if not diff.empty:
            out.write("Differences found:\n")
            out.write(diff.to_string(index=False))
        else:
            out.write("The files match perfectly.\n")


file1 = 'egm_GCR_Moritz.txt'


compare_columns(file1)

