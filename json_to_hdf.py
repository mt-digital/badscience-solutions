import sys

from util import ExperimentData


def main(jsons_dir, dest_hdf_name):

    ExperimentData.from_jsons(jsons_dir, dest_hdf_name)


if __name__ == "__main__":

    if sys.argv[1] == "-h":
        print("python json_to_hdf.py <jsons_dir> <dest_hdf_name>")
        sys.exit(0)

    main(sys.argv[1], sys.argv[2])
