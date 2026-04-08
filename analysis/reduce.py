import os
from coffea.util import load, save
from helpers.futures_patch import patch_mp_connection_bpo_17560


def _parse_csv_set(s):
    if s is None:
        return None
    vals = [x.strip() for x in s.split(",") if x.strip()]
    return set(vals) if vals else None


def _match_dataset(name, dataset_set):
    if dataset_set is None:
        return True
    return any(d in name for d in dataset_set)


def _match_exclude(name, exclude_set):
    if exclude_set is None:
        return False
    return any(d in name for d in exclude_set)


def reduce(folder, _dataset=None, _exclude=None, variable=None):
    dataset_set = _parse_csv_set(_dataset)
    exclude_set = _parse_csv_set(_exclude)
    variable_set = _parse_csv_set(variable)

    # group .futures files by pdi = filename.split("____")[0]
    groups = {}
    with os.scandir(folder) as it:
        for entry in it:
            if not entry.is_file():
                continue
            if ".futures" not in entry.name:
                continue
            pdi = entry.name.split("____")[0]
            groups.setdefault(pdi, []).append(entry.path)

    for pdi, files in groups.items():
        if not _match_dataset(pdi, dataset_set):
            continue
        if _match_exclude(pdi, exclude_set):
            continue

        print(f"[dataset] {pdi}  ({len(files)} files)")

        acc = {}

        for i, filename in enumerate(files, 1):
            print(f"  [{i}/{len(files)}] loading {os.path.basename(filename)}")
            hin = load(filename)

            for k, h in hin.items():
                if variable_set is not None and k not in variable_set:
                    continue

                if k not in acc:
                    acc[k] = h
                else:
                    acc[k] = acc[k] + h

            del hin

        for k, h in acc.items():
            out = {k: h}
            outname = os.path.join(folder, f"{k}--{pdi}.reduced")
            print(f"  saving {outname}")
            save(out, outname)


if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-f", "--folder", help="folder", dest="folder")
    parser.add_option("-d", "--dataset", help="dataset", dest="dataset", default=None)
    parser.add_option("-e", "--exclude", help="exclude", dest="exclude", default=None)
    parser.add_option("-v", "--variable", help="variable", dest="variable", default=None)
    (options, args) = parser.parse_args()

    patch_mp_connection_bpo_17560()
    reduce(options.folder, options.dataset, options.exclude, options.variable)