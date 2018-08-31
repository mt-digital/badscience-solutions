import numpy as np
import json
import os

from glob import glob
from h5py import File


class ExperimentData:

    def __init__(self, hdf_path):
        hdf = File(hdf_path, "r")
        self.policies = list(hdf.keys())
        self.award_amounts = list(
            hdf[self.policies[0]].keys()
        )
        self.negres_rates = list(
            hdf[self.policies[0]][self.award_amounts[0]].keys()
        )
        self.fpdrs = list(
            hdf[self.policies[0]][
                    self.award_amounts[0]][
                        self.negres_rates[0]].keys()
        )

        self.hdf = hdf

    @classmethod
    def from_jsons(cls, jsons_dir, hdf_path):
        """
        Create new instance from a directory of jsons (raw output from
        scimod-agency). Write the new HDF to hdf_path. HDF stays open unless
        closed explicitly. Using Python destructor to close in case user
        forgets.
        """
        hdf = File(hdf_path, "w")

        json_files = glob(os.path.join(jsons_dir, "*"))
        for jf in json_files:
            j = json.load(open(jf, "r"))
            md = j["metadata"]
            policy = md["policy"]

            parameters = md["parameters"]
            award_amount = parameters["awardAmount"]
            pubneg_rate = parameters["publishNegativeResultRate"]
            fpdr = parameters["falsePositiveDetectionRate"]

            group = "/{}/{}/{:.2f}/{:.2f}".format(
                policy, award_amount, pubneg_rate, fpdr
            )
            try:
                hdf.create_group(group)
                for measure in [
                            'falseDiscoveryRate',
                            'falsePositiveRate',
                            'nPublications'
                        ]:
                    try:
                        hdf[group].create_dataset(
                            measure, data=np.array(j[measure], dtype='float'),
                            compression="gzip", compression_opts=9)
                    except KeyError:
                        pass

            except ValueError:
                pass

        hdf.close()

        return cls(hdf_path)

    def __del__(self):
        self.hdf.close()
