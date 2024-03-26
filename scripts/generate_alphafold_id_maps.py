from pathlib import Path
from tqdm import tqdm

if __name__ == "__main__":
    afdb_path = Path("/vol/structure/alphafold/taxa")
    output_map = Path("/home/khippe/afdb_id_map.tsv")

    with open(output_map, "w") as f:
        f.write("id\tpath\n")
        for fpath in tqdm(afdb_path.glob("**/*.*")):
            if fpath.is_file() and "model" in fpath.stem:
                f.write(f"{fpath.stem}\t{fpath}\n")
